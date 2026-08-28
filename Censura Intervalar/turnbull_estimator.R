### ==========================================================================
### ESTIMADOR DE TURNBULL (NPMLE) PARA DADOS COM CENSURA INTERVALAR
### ==========================================================================
#
# Baseado na verossimilhança geral da sua anotação:
#
#     L(theta) = prod_{i=1}^{n} [ S(L_i) - S(U_i) ]
#
# Essa forma comporta todos os tipos de censura que você listou:
#   T = t          (exato)            -> L_i = U_i = t   (use um epsilon pequeno, ver nota)
#   T in [L_i,U_i] (intervalar)       -> L_i, U_i finitos
#   T > L_i        (censura à direita)-> U_i = Inf
#   T < U_i        (censura à esquerda)-> L_i = 0
#
# Como S(.) é desconhecida e não-paramétrica, Turnbull (1976) mostrou que o
# MLE de S só pode colocar massa em um conjunto finito de intervalos
# "não decomponíveis" (os chamados intervalos de Turnbull), e que o problema
# se resolve por um algoritmo EM / autoconsistente. É exatamente o que as
# funções abaixo fazem.
#

## ----------------------------------------------------------------------
## PARTE 1 — Encontrar os intervalos de Turnbull (q_j, p_j]
## ----------------------------------------------------------------------
## Ideia: junta todos os L's e U's, ordena, e forma intervalos consecutivos.
## Só ficam os intervalos que estão *contidos* em pelo menos uma observação
## (L_i, U_i] — são esses os únicos lugares onde a verossimilhança "enxerga"
## massa de probabilidade (equivalente a S(L_i) - S(U_i) = soma das massas
## s_j dos intervalos de Turnbull contidos em (L_i,U_i]).

find_turnbull_intervals <- function(L, U) {
  pts <- sort(unique(c(L, U)))
  m <- length(pts)
  q <- pts[-m]          # extremos esquerdos candidatos
  p <- pts[-1]          # extremos direitos candidatos
  keep <- sapply(seq_along(q), function(j) any(L <= q[j] & U >= p[j]))
  data.frame(q = q[keep], p = p[keep])
}


## ----------------------------------------------------------------------
## PARTE 2 — Algoritmo EM / autoconsistente (o coração do estimador)
## ----------------------------------------------------------------------
## alpha[i,j] = 1{ (q_j,p_j] esta contido em (L_i,U_i] }
##   -> é o "indicador" que aparece implicitamente na sua fórmula:
##      S(L_i) - S(U_i) = sum_j alpha[i,j] * s_j
##
## Passo E: distribui cada observação i entre os intervalos de Turnbull
##          proporcionalmente à massa atual s_j (só entre os intervalos
##          compatíveis, alpha[i,j]=1)
## Passo M: atualiza s_j = média das "responsabilidades" da coluna j
##
## Isso maximiza exatamente L(theta) = prod_i [S(L_i)-S(U_i)] = prod_i sum_j alpha_ij s_j

turnbull_em <- function(L, U, tol = 1e-9, maxit = 10000) {
  # Tratamento para tempos exatos (L_i == U_i)
  exact_idx <- (L == U)
  if (any(exact_idx)) {
    U[exact_idx] <- L[exact_idx] + 1e-7
  }
  ints <- find_turnbull_intervals(L, U)
  q <- ints$q; p <- ints$p
  k <- length(q); n <- length(L)

  alpha <- matrix(0, n, k)
  for (i in 1:n) alpha[i, ] <- as.numeric(L[i] <= q & U[i] >= p)

  s <- rep(1 / k, k)   # chute inicial: massa uniforme
  for (iter in 1:maxit) {
    numer <- sweep(alpha, 2, s, `*`)          # alpha_ij * s_j
    denom <- rowSums(numer)                    # S(L_i)-S(U_i) da iteração atual
    denom[denom == 0] <- 1e-12
    mu <- numer / denom                         # passo E (responsabilidades)
    mu[is.na(mu)] <- 0
    # Normalização preventiva de contingência
    s_new <- colSums(mu) / n
    total_s <- sum(s_new)
    if (total_s > 0) s_new <- s_new / total_s
    
    # Checagem de convergência segura contra NAs
    diff <- abs(s_new - s)
    if (all(!is.na(diff)) && max(diff) < tol) { 
      s <- s_new
      break 
    }
    s <- s_new
  }
  list(q = q, p = p, s = s, iter = iter, mu = mu, alpha = alpha, n = n)
}


## ----------------------------------------------------------------------
## PARTE 3 — Tabela no estilo do summary(survfit(...)):
##            times, n.risk, n.event, survival, std.err, lower/upper 95% CI
## ----------------------------------------------------------------------
## Turnbull mostrou que a solução autoconsistente coincide com um
## produto-limite "generalizado":
##     d_j (n.event) = n * s_j                (nº esperado de eventos no intervalo j)
##     n_j (n.risk)  = soma de d_k para k >= j  (nº esperado ainda "em risco")
##     S(p_j) = prod_{k<=j} (1 - d_k/n_k)      <=> 1 - cumsum(s)   (mesma coisa)
## O erro-padrão usa a fórmula de Greenwood aplicada a esses d_j, n_j
## "pseudo" (é a mesma lógica do Kaplan-Meier, só que sobre os intervalos
## de Turnbull ao invés de tempos exatos).

turnbull_table <- function(fit, conf.level = 0.95) {
  n <- fit$n
  d <- n * fit$s
  n_risk <- rev(cumsum(rev(d)))
  surv <- 1 - cumsum(fit$s)

  # Greenwood
  term <- ifelse(n_risk > d, d / (n_risk * (n_risk - d)), 0)
  cum_var <- cumsum(term)
  se <- surv * sqrt(cum_var)

  # IC por transformação log-log (mesma default do survival::survfit),
  # evita limites fora de [0,1]
  z <- qnorm(1 - (1 - conf.level) / 2)
  theta <- exp(z * se / (surv * log(surv)))
  theta[!is.finite(theta)] <- 1
  lower <- surv^(1 / theta)
  upper <- surv^theta
  edge <- surv %in% c(0, 1)
  lower[edge] <- surv[edge]; upper[edge] <- surv[edge]

  data.frame(
    q = fit$q, p = fit$p,
    n.risk   = round(n_risk, 3),
    n.event  = round(d, 3),
    survival = round(surv, 4),
    std.err  = round(se, 4),
    lower95  = round(pmax(lower, 0), 4),
    upper95  = round(pmin(upper, 1), 4)
  )
}


## ----------------------------------------------------------------------
## PARTE 3b (opcional) — IC por bootstrap
## ----------------------------------------------------------------------
## A fórmula de Greenwood acima é uma aproximação (nem sempre boa com poucos
## dados). O jeito "padrão-ouro" na literatura de Turnbull é bootstrap
## reamostrando os indivíduos e refazendo o EM.

turnbull_boot_ci <- function(L, U, B = 500, conf.level = 0.95, times = NULL) {
  n <- length(L)
  fit0 <- turnbull_em(L, U)
  if (is.null(times)) times <- fit0$p

  boot_surv <- matrix(NA, B, length(times))
  for (b in 1:B) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    fb <- turnbull_em(L[idx], U[idx])
    Sb <- 1 - cumsum(fb$s)
    boot_surv[b, ] <- stepfun(fb$p, c(1, Sb))(times)
  }
  alpha <- 1 - conf.level
  data.frame(
    time  = times,
    lower = apply(boot_surv, 2, quantile, probs = alpha / 2, na.rm = TRUE),
    upper = apply(boot_surv, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
}


### NOTA — encaixando o modelo S(t) = S0(t)^exp(x'beta) (Cox) da sua folha
### ==========================================================================
## O que está acima é a parte NÃO-PARAMÉTRICA (estima S0 livre, sem covariável).
## Para colocar covariáveis x'beta como na sua última linha, o caminho natural
## é usar a MESMA log-verossimilhança de intervalo, mas trocando S(L_i), S(U_i)
## por S0(L_i)^exp(x_i'beta), S0(U_i)^exp(x_i'beta) e maximizando com relação a
## (S0 nos intervalos de Turnbull, beta) simultaneamente — é o chamado modelo
## de Cox com censura intervalar (ex.: pacote `icenReg`, função `ic_sp`).
## Se for esse o próximo passo do seu trabalho, me avisa que eu te mostro como
## estender o EM acima para incluir o beta (via Newton-Raphson dentro do M-step).
