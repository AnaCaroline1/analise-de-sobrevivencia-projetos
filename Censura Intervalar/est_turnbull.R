turnbull <- function(L, U, group = NULL, conf.level = 0.95,
                     maxit = 2000, tol = 1e-9) {
  
  stopifnot(length(L) == length(U))
  
  if (is.null(group)) {
    group <- factor(rep("Overall", length(L)))
  } else {
    group <- factor(group)
  }
  
  fits <- lapply(levels(group), function(g) {
    idx <- group == g
    .turnbull_fit(L[idx], U[idx], conf.level = conf.level,
                  maxit = maxit, tol = tol)
  })
  names(fits) <- levels(group)
  
  lbl <- paste0(conf.level, "LCL")
  ubl <- paste0(conf.level, "UCL")
  
  tab <- do.call(rbind, lapply(names(fits), function(g) {
    f <- fits[[g]]
    data.frame(grupo = g, records = f$n, n.events = f$n.events,
               median = f$median, LCL = f$lcl, UCL = f$ucl)
  }))
  names(tab)[names(tab) == "LCL"] <- lbl
  names(tab)[names(tab) == "UCL"] <- ubl
  if (nlevels(group) == 1) tab$grupo <- NULL
  
  structure(list(table = tab, fits = fits, conf.level = conf.level),
            class = "turnbull")
}

print.turnbull <- function(x, digits = 3, ...) {
  print(x$table, row.names = FALSE, digits = digits)
  invisible(x)
}

# curva S(t) estimada, empilhada por grupo (para plot/inspecao)
survival_table <- function(x, ...) UseMethod("survival_table")
survival_table.turnbull <- function(x, ...) {
  do.call(rbind, lapply(names(x$fits), function(g) {
    d <- x$fits[[g]]$table_S
    d$grupo <- g
    d
  }))
}

# -------------------------------------------------------------
# Ajuste interno para um unico grupo
# -------------------------------------------------------------
.turnbull_fit <- function(L, U, conf.level = 0.95, maxit = 2000, tol = 1e-9) {
  
  n  <- length(L)
  Ui <- ifelse(is.na(U) | is.infinite(U), Inf, U)
  delta <- as.integer(is.finite(Ui))   # 1 = intervalo finito, 0 = censura a direita
  
  # pontos de tempo candidatos: uniao dos limites observados (L e U finitos)
  tau <- sort(unique(c(0, L, Ui[is.finite(Ui)])))
  m   <- length(tau) - 1               # numero de intervalos (tau_{j-1}, tau_j]
  
  lo <- tau[-length(tau)]              # tau_{j-1}, vetor de tamanho m
  hi <- tau[-1]                        # tau_j
  
  # -- matriz de indicadores alpha_ij --
  # Regra geral: (tau_{j-1}, tau_j] contido em (L_i, U_i]?
  # CASO ESPECIAL -- tempo exato de falha (L_i = U_i = t): o intervalo
  # (L_i, U_i] = (t, t] e VAZIO, entao a regra geral zeraria alpha_ij para
  # todo j e a observacao desapareceria da verossimilhanca (perda de
  # massa). Isso corresponde a linha "T = t -> f(t)" da tabela: aqui o
  # individuo so pode ter "causado" o evento no unico tau_j que coincide
  # com o proprio tempo t, entao tratamos esse caso separadamente.
  exato <- is.finite(Ui) & (Ui == L)
  
  alpha <- matrix(0L, n, m)
  for (i in seq_len(n)) {
    if (exato[i]) {
      alpha[i, ] <- as.integer(hi == L[i])
    } else {
      alpha[i, ] <- as.integer(lo >= L[i] & hi <= Ui[i])
    }
  }
  
  # valor inicial para S: decaimento linear entre 1 e 0
  S <- seq(1, 0, length.out = m + 1)
  
  d <- rep(0, m); n_risco <- rep(0, m)
  
  for (it in seq_len(maxit)) {
    
    p <- -diff(S)
    p[p < 0] <- 0
    
    peso <- as.vector(alpha %*% p)          # soma_k alpha_ik p_k, por individuo
    peso[peso == 0] <- Inf                  # evita divisao por zero
    
    num <- alpha * matrix(p, n, m, byrow = TRUE)
    d   <- colSums(num / matrix(peso, n, m))
    
    n_risco <- rev(cumsum(rev(d)))          # n_j = soma_{k=j}^{m} d_k
    fator   <- ifelse(n_risco > 0, 1 - d / n_risco, 1)
    S_novo  <- c(1, cumprod(fator))
    
    if (max(abs(S_novo - S)) < tol) { S <- S_novo; break }
    S <- S_novo
  }
  
  # log-verossimilhanca (para conferencia / diagnostico)
  Sfun <- function(t) {
    tt <- ifelse(is.infinite(t), max(tau), t)
    s  <- approx(tau, S, xout = tt, method = "constant", rule = 2, f = 0)$y
    ifelse(is.infinite(t), 0, s)
  }
  ll <- sum(delta * log(pmax(Sfun(L) - Sfun(Ui), 1e-12)) +
              (1 - delta) * log(pmax(Sfun(L), 1e-12)))
  
  # erro-padrao de log(-log S) via formula de Greenwood (transformacao log-log)
  var_term <- ifelse(n_risco > d & n_risco > 0,
                     d / (n_risco * (n_risco - d)), 0)
  Spos <- pmin(pmax(S[-1], 1e-10), 1 - 1e-10)
  var_cloglog <- cumsum(var_term) / (log(Spos))^2
  se_cloglog  <- sqrt(pmax(var_cloglog, 0))
  
  z   <- qnorm(1 - (1 - conf.level) / 2)
  cll <- log(-log(Spos))
  lower <- exp(-exp(cll + z * se_cloglog))   # banda inferior de S(t)
  upper <- exp(-exp(cll - z * se_cloglog))   # banda superior de S(t)
  
  # mediana e seu IC: tempo em que a curva de S(t) cruza 0.5
  cruza <- function(tempos, surv) {
    idx <- which(surv <= 0.5)
    if (length(idx) == 0) return(NA_real_)
    tempos[min(idx)]
  }
  mediana <- cruza(hi, S[-1])
  lcl_med <- cruza(hi, lower)   # curva inferior cruza 0.5 mais cedo -> limite inferior da mediana
  ucl_med <- cruza(hi, upper)   # curva superior cruza 0.5 mais tarde -> limite superior da mediana
  
  list(
    tau = tau, S = S, d = d, n_risco = n_risco,
    n = n, n.events = sum(delta), loglik = ll,
    median = mediana, lcl = lcl_med, ucl = ucl_med,
    table_S = data.frame(time = hi, S = S[-1], lower = lower, upper = upper)
  )
}