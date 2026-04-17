log_vero_surv_weibull_R <- function(theta, t, status) {
  
  # theta[1] = alpha (scale)
  # theta[2] = gamma (shape)
  alpha <- theta[1]
  gamma <- theta[2]
  
  # Evita valores não positivos para os parâmetros durante a otimização
  if (alpha <= 0 || gamma <= 0) return(1e10)
  
  # log(f(t)): log-densidade para os eventos (status == 1)
  log_f <- dweibull(t, shape = gamma, scale = alpha, log = TRUE)
  
  # log(S(t)): log-sobrevivência para as censuras (status == 0)
  # lower.tail = FALSE retorna P(X > x), e log.p = TRUE retorna o log disso
  log_S <- pweibull(t, shape = gamma, scale = alpha, lower.tail = FALSE, log.p = TRUE)
  
  # Soma ponderada pelo status:
  # Se status=1, contribui log_f. Se status=0, contribui log_S.
  Log_L <- sum(status * log_f + (1 - status) * log_S)
  
  if (!is.finite(Log_L)) return(1e10)
  
  return(-Log_L) # Retorna negativo para minimizar
}

tempos  
cens

ajust3 <- optim(par = c(1, 1), log_vero_surv_weibull_R,
                t = tempos, status = cens,
                method = "BFGS",hessian = T)
#############################################

ekm <- function(tempos,status){
  ordem <- order(tempos)
  t <- tempos[ordem]
  s <- status[ordem]
  
  n <- length(t)
  surv <- numeric(n)
  p_surv <- 1
  
  for(i in 1:n) {
    # n_risco: quantos restam (incluindo o atual)
    n_risco <- n - i + 1
    # Se houve evento, a probabilidade condicional de sobreviver é (n-1)/n
    if(s[i] == 1) {
      p_surv <- p_surv * ((n_risco - 1) / n_risco)
    }
    surv[i] <- p_surv
  }
  
  return(list(tempo = t, sobrevivencia = surv))
}
km_c_bx <- ekm(tempos,cens)
km_c_bx

##########################################################

# 1. Não filtre os valores iguais a 1. 
# Se você precisa filtrar apenas para evitar S(t)=0 (por causa de logs, por exemplo), faça:
indices <- which(km_c_bx$sobrevivencia > 0) 
t_km <- km_c_bx$tempo[indices]
s_km <- km_c_bx$sobrevivencia[indices]

# 2. Se o seu km_c_bx não tiver o tempo 0, adicione-o manualmente para o plot:
t_plot <- c(0, t_km)
s_plot <- c(1, s_km)

# 3. Parametros e t_seq (mantenha como está, t_seq começando em 0 é ótimo)
alpha_est <- ajust3$par[1]
gamma_est <- ajust3$par[2]
t_seq <- seq(0, max(t_plot), length.out = 100)
s_weibull <- pweibull(t_seq, shape = gamma_est, scale = alpha_est, lower.tail = FALSE)

# 4. Plotar usando os vetores com o tempo zero incluído
plot(t_plot, s_plot, type = "s", lwd = 2,
       xlab = "Tempo (t)", ylab = "S(t)", 
       main = "Ajuste Weibull vs. Kaplan-Meier",
       xlim = c(0, max(t_plot)), 
       ylim = c(0, 1),
       xaxs = "i", yaxs = "i") # 'i' faz os eixos encostarem exatamente no 0 e 1
  
  lines(t_seq, s_weibull, col = "darkred", lwd = 2, lty = 2)
  
  legend("topright", legend = c("Kaplan-Meier", "Weibull Ajustada"),
         col = c("black", "darkred"), lty = c(1, 2), lwd = 2)

# 2. Transformações lineares
y_weibull <- log(-log(s_km))
x_weibull <- log(t_km)

# 3. Plotar os dados observados
plot(x_weibull, y_weibull, pch = 16, col = "purple",
     xlab = "log(t)", ylab = "log(-log(S(t)))",
     main = "Gráfico de Linearização Weibull")

# 4. Adicionar a reta teórica baseada no seu ajuste
# Intercepto: -gamma * log(alpha) | Inclinação: gamma
abline(a = -gamma_est * log(alpha_est), b = gamma_est, col = "red", lwd = 2)


# Calculando os intervalos de confiança

Sigma <- solve(ajust3$hessian) 

# 1. Definir tempos para o plot
t_seq <- seq(0, max(tempos), length.out = 100)

# 2. Calcular S(t) e o Erro Padrão para cada t
s_pred <- numeric(length(t_seq))
se_s   <- numeric(length(t_seq))

for(i in 1:length(t_seq)) {
  ti <- t_seq[i]
  # S(t) estimado
  s_pred[i] <- exp(-(ti/alpha_est)^gamma_est)
  
  # Gradiente no ponto ti
  d_alpha <- s_pred[i] * gamma_est * (ti^gamma_est / alpha_est^(gamma_est+1))
  d_gamma <- -s_pred[i] * (ti/alpha_est)^gamma_est * log(ti/alpha_est + 1e-6)
  grad <- c(d_alpha, d_gamma)
  
  # Variância via Delta Method
  var_s <- t(grad) %*% Sigma %*% grad
  se_s[i] <- sqrt(var_s)
}

# 3. Limites do IC (95%)
ic_inf <- pmax(0, s_pred - 1.96 * se_s) # Garante que não seja < 0
ic_sup <- pmin(1, s_pred + 1.96 * se_s) # Garante que não seja > 1

plot(t_plot, s_plot, type = "s", lwd = 2,
     xlab = "Tempo (t)", ylab = "S(t)", 
     main = "Ajuste Weibull vs. Kaplan-Meier",
     xlim = c(0, max(t_plot)), 
     ylim = c(0, 1),
     xaxs = "i",
     yaxs = "i"
     ) # 'i' faz os eixos encostarem exatamente no 0 e 1

lines(t_seq, s_weibull, col = "darkred", lwd = 2, lty = 2)
lines(t_seq, ic_inf, col = "red", lty = 2)
lines(t_seq, ic_sup, col = "red", lty = 2)

# Opcional: Polígono de preenchimento
polygon(c(t_seq, rev(t_seq)), c(ic_inf, rev(ic_sup)), 
        col = rgb(1, 0, 0, 0.2), border = NA)
legend("topright", legend = c("Kaplan-Meier", "Weibull Ajustada"),
       col = c("black", "darkred"), lty = c(1, 2), lwd = 2)
