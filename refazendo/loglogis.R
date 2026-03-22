S_loglogistica <- function(theta, t, delta) {
  
  alpha <- exp(theta[1])
  gamma <- exp(theta[2])
  
  log_f <- log(gamma) - (gamma * log(alpha)) + (gamma - 1) * log(t) - 2 * log(1 + (t/alpha)^gamma)
  
  log_S <- -log(1 + (t/alpha)^gamma)
  
  LOG_L <- sum(delta * log_f + (1 - delta) * log_S)
  
  return(-LOG_L)
}

inicio <- c(log(mean(tempos)), 0)

fit_ll <- optim(par = inicio, fn = S_loglogistica, t = tempos, 
                delta = cens, method = "BFGS")
fit_ll

km_c_bx
# 1. Extração correta dos parâmetros (voltando para a escala original)
alpha_est <- exp(fit_ll$par[1])
gamma_est <- exp(fit_ll$par[2])

# 2. Sequência de tempo começando do zero para o gráfico encostar no eixo Y
t_seq <- seq(0, max(tempos, na.rm = TRUE), length.out = 100)

# 3. Calcular a sobrevivência teórica CORRETA
s_loglogis <- 1 / (1 + (t_seq / alpha_est)^gamma_est)

# 4. Plotagem com os ajustes de limites que discutimos antes
plot(km_c_bx$tempo, km_c_bx$sobrevivencia, 
     type = "s", lwd = 2,
     xlab = "Tempo (t)", ylab = "S(t)", 
     main = "Ajuste Log-Logístico vs. Kaplan-Meier",
     xlim = c(0, max(tempos)), # Força começar no t=0
     ylim = c(0, 1),
     xaxs = "i", yaxs = "i")   # Estética profissional (encosta nos eixos)

# Adiciona a curva teórica
lines(t_seq, s_loglogis, col = "darkgreen", lwd = 2, lty = 2)

# Adiciona a legenda
legend("topright", 
       legend = c("Kaplan-Meier", "Log-Logística Ajustada"),
       col = c("black", "darkgreen"), 
       lty = c(1, 2), 
       lwd = 2,
       bty = "n")

# 1. Preparar os dados do Kaplan-Meier
# Precisamos filtrar S(t) < 1 e S(t) > 0 para o logito não dar infinito
indices_lin <- which(km_c_bx$sobrevivencia > 0 & km_c_bx$sobrevivencia < 1)
t_km <- km_c_bx$tempo[indices_lin]
s_km <- km_c_bx$sobrevivencia[indices_lin]

# 2. Definir as coordenadas de linearização
y_lin <- log((1 - s_km) / s_km)
x_lin <- log(t_km)

# 3. Criar o gráfico de dispersão
plot(x_lin, y_lin, 
     pch = 16, col = "darkblue",
     xlab = "log(t)", 
     ylab = "log((1 - S(t)) / S(t))",
     main = "Gráfico de Linearização: Log-Logística")

# 4. Adicionar a reta teórica baseada nos seus parâmetros estimados (fit_ll)
# Lembre-se: alpha_est e gamma_est devem ser os valores na escala original
alpha_est <- exp(fit_ll$par[1])
gamma_est <- exp(fit_ll$par[2])

# O intercepto é -gamma * log(alpha) e a inclinação é gamma
abline(a = -gamma_est * log(alpha_est), b = gamma_est, 
       col = "red", lwd = 2, lty = 2)

grid()