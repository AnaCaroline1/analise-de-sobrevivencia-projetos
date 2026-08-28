# =========================================================================
# ESTIMADOR DE TURNBULL (NPMLE) PARA CENSURA INTERVALAR
# Com extensao para grupos, tabela estilo survfit() e IC via bootstrap
# =========================================================================
#
# Liga direto com a foto:
#   T ~ F(.|theta),  T in [L,U]
#   L(theta) = prod_{i=1}^{n} [S(Li) - S(Ui)]^{delta_i}
#
# Como S(.) e' modelada de forma nao-parametrica (massas p_k em subintervalos
# de Turnbull), temos S(Li)-S(Ui) = soma das massas p_k contidas em [Li,Ui]
# = (A %*% p)_i.  Ou seja, a verossimilhanca vira:
#
#   L(p) = prod_i (A p)_i        sujeito a  sum(p) = 1, p >= 0
#
# O algoritmo de Turnbull (self-consistency) e' um EM que maximiza essa L(p).
# =========================================================================

suppressMessages(library(survival))

# --------------------------------------------------------------------
# 1) GRADE DE TEMPOS (tau) -- pontos candidatos para massa de probabilidade
# --------------------------------------------------------------------
cria.tau <- function(data){
  l <- data$left
  r <- data$right
  tau <- sort(unique(c(l, r[is.finite(r)])))
  return(tau)
}

# --------------------------------------------------------------------
# 2) CHUTE INICIAL de p (so' para dar ponto de partida ao EM)
# --------------------------------------------------------------------
S.ini <- function(tau){
  m <- length(tau)
  ekm <- survfit(Surv(tau[1:m-1], rep(1, m-1)) ~ 1)
  So <- c(1, ekm$surv)
  p <- -diff(So)
  return(p)
}

# --------------------------------------------------------------------
# 3) MATRIZ A -- liga cada individuo aos subintervalos (tau_{k-1},tau_k]
#    compativeis com [Li,Ui]. A%*%p = S(Li)-S(Ui) = termo da verossimilhanca
# --------------------------------------------------------------------
cria.A <- function(data, tau){
  tau12 <- cbind(tau[-length(tau)], tau[-1])
  interv <- function(x, inf, sup) ifelse(x[1] >= inf & x[2] <= sup, 1, 0)
  A <- apply(tau12, 1, interv, inf = data$left, sup = data$right)
  id.lin.zero <- which(apply(A == 0, 1, all))
  if (length(id.lin.zero) > 0) A <- A[-id.lin.zero, ]
  return(A)
}

# --------------------------------------------------------------------
# 4) ALGORITMO EM (self-consistency) -- maximiza log L(p) = sum(log(A p))
#    Passo E: C_i = (Ap)_i = S(Li)-S(Ui)   (contribuicao de cada individuo)
#    Passo M: p_k <- p_k * sum_i[A_ik / C_i] / n
# --------------------------------------------------------------------
turnbull_em <- function(p, A, eps = 1e-6, iter.max = 2000){
  n <- nrow(A)
  Q <- matrix(1, ncol(A))
  iter <- 0
  repeat{
    iter <- iter + 1
    Q <- p
    C <- A %*% p                        # C_i = S(Li)-S(Ui)  -> verossimilhanca individual
    C <- pmax(C, 1e-10)                 # PISO NUMERICO: evita divisao por zero / NaN
                                         # quando algum individuo fica com C_i=0
    p <- p * ((t(A) %*% (1 / C)) / n)   # passo M do EM
    maxdiff <- max(abs(p - Q))
    if (!is.finite(maxdiff)) {
      stop("EM instavel (maxdiff nao finito) na iteracao ", iter,
           " - dados degenerados (poucos subintervalos validos).")
    }
    if (maxdiff < eps | iter >= iter.max) break
  }
  loglik <- sum(log(pmax(A %*% p, 1e-10)))
  list(p = as.vector(p), iter = iter, maxdiff = maxdiff, loglik = loglik)
}

# --------------------------------------------------------------------
# 4b) PRE-PROCESSAMENTO: trata tempos exatos (left==right)
#    O algoritmo basico de Turnbull so lida com censura intervalar
#    "verdadeira" (L<U). Um tempo exato T=t precisa virar uma janela
#    bem estreita (t-eps, t] para caber na grade de subintervalos.
# --------------------------------------------------------------------
ajustar_exatos <- function(data, eps = 1e-3){
  exatos <- which(data$left == data$right)
  if (length(exatos) > 0){
    data$left[exatos] <- data$left[exatos] - eps
  }
  data
}

# --------------------------------------------------------------------
# 5) AJUSTE COMPLETO -- roda 1-4 e monta tabela estilo survfit()
#    n.event_k  = massa esperada de eventos no intervalo k  ~ round(p_k * n)
#    n.risk_k   = individuos "em risco" antes do intervalo k
#    (convencao tipo tabua de vida -- nao ha risk-set exato p/ dado intervalar)
# --------------------------------------------------------------------
turnbull_fit <- function(data, eps = 1e-6, iter.max = 2000, init = c("uniform","km")){
  init <- match.arg(init)
  data <- ajustar_exatos(data)
  n <- nrow(data)
  tau <- cria.tau(data)
  A <- cria.A(data, tau)                # nrow(A)=individuos (validos), ncol(A)=length(tau)-1
  if (is.null(dim(A)) || ncol(A) < 1){
    stop("Dados degenerados: nenhum subintervalo valido (poucos valores distintos ",
         "de left/right nesta amostra).")
  }
  tau_ef <- tau[-1]                     # limite superior de cada subintervalo (tau_{k-1},tau_k]

  # Inicializacao UNIFORME por padrao: como o EM e' multiplicativo, comecar
  # com p_k=0 em algum lugar (o que S.ini() pode fazer) trava esse componente
  # em zero para sempre. Uniforme evita esse problema.
  p0 <- if (init == "uniform") rep(1 / ncol(A), ncol(A)) else S.ini(tau)

  fit <- turnbull_em(p0, A, eps = eps, iter.max = iter.max)
  p <- fit$p

  surv <- c(1, 1 - cumsum(p))
  n.event <- round(p * n)
  n.risk  <- n - c(0, cumsum(n.event))[-(length(n.event) + 1)]

  out <- data.frame(
    time    = tau_ef,
    n.risk  = n.risk,
    n.event = n.event,
    p.mass  = p,
    surv    = surv[-1]
  )
  attr(out, "loglik") <- fit$loglik
  attr(out, "iter")   <- fit$iter
  attr(out, "n")      <- n
  out
}

# --------------------------------------------------------------------
# 6) STD.ERR e IC 95% via BOOTSTRAP NAO-PARAMETRICO
#    (reamostra individuos, reajusta Turnbull, guarda S(t) em cada tempo
#     da grade original -- funcao "step" para casar tempos entre replicas)
# --------------------------------------------------------------------
turnbull_boot <- function(data, B = 300, seed = 1){
  set.seed(seed)
  base <- turnbull_fit(data)   # se os dados originais forem degenerados, erro aqui e' esperado/util
  times <- base$time
  n <- nrow(data)

  surv_at <- function(fit_i, t){
    if (length(fit_i$time) < 1) return(rep(NA_real_, length(t)))
    # funcao S(t) em degrau a partir do ajuste bootstrap fit_i
    sfun <- stepfun(fit_i$time, c(1, fit_i$surv))
    sfun(t)
  }

  S_boot <- matrix(NA, nrow = B, ncol = length(times))
  b <- 1
  tent <- 0
  falhas <- 0
  while (b <= B && tent < B * 5){
    tent <- tent + 1
    idx <- sample(seq_len(n), n, replace = TRUE)
    d_b <- data[idx, ]
    fit_b <- tryCatch(turnbull_fit(d_b), error = function(e) NULL)
    if (!is.null(fit_b) && length(fit_b$time) >= 1){
      linha <- tryCatch(surv_at(fit_b, times), error = function(e) NULL)
      if (!is.null(linha)){
        S_boot[b, ] <- linha
        b <- b + 1
        next
      }
    }
    falhas <- falhas + 1
  }
  if (falhas > 0){
    message(falhas, " de ", tent, " reamostragens bootstrap foram descartadas ",
            "(dados degenerados na reamostra) - normal em amostras pequenas.")
  }
  S_boot <- S_boot[seq_len(b - 1), , drop = FALSE]
  if (nrow(S_boot) < 2){
    stop("Bootstrap falhou na quase totalidade das reamostragens - amostra ",
         "provavelmente pequena/degenerada demais para este B. Tente reduzir B ",
         "ou aumentar o eps do EM.")
  }

  se <- apply(S_boot, 2, sd, na.rm = TRUE)

  # IC 95% na escala log(-log(S)) (mesma ideia usada no Kaplan-Meier),
  # que mantem o intervalo dentro de [0,1]
  S <- base$surv
  z <- qnorm(0.975)
  theta <- log(-log(S))
  se_theta <- se / (S * abs(log(S)))
  se_theta[!is.finite(se_theta)] <- NA

  lower <- exp(-exp(theta + z * se_theta))
  upper <- exp(-exp(theta - z * se_theta))

  base$std.err <- se
  base$lower95 <- pmin(lower, 1, na.rm = FALSE)
  base$upper95 <- pmin(upper, 1, na.rm = FALSE)
  base
}

# --------------------------------------------------------------------
# 7) EXTENSAO PARA GRUPOS -- roda tudo (4-6) separadamente por grupo
#    (equivalente a survfit(Surv(...) ~ grupo, type="interval2"))
# --------------------------------------------------------------------
turnbull_by_group <- function(data, group_col = "grupo", B = 300, seed = 1){
  if (!group_col %in% names(data)){
    stop("Coluna '", group_col, "' nao existe em 'data'. Colunas disponiveis: ",
         paste(names(data), collapse = ", "),
         ". Use o argumento group_col para indicar o nome correto.")
  }
  grupos <- split(data, data[[group_col]])
  res <- lapply(names(grupos), function(g){
    tab <- turnbull_boot(grupos[[g]], B = B, seed = seed)
    tab$grupo <- g
    tab
  })
  do.call(rbind, res)
}
