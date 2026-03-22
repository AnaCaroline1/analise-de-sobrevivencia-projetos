log_vero_surv_lnorm <- function(theta,t,status){
  mu    <- theta[1]
  sigma <- theta[2]
  log_S <- plnorm(t, meanlog = mu, sdlog = sigma, lower.tail = FALSE, log.p = TRUE)
  
  log_f <- -0.5*log(2*pi) - log(t) - log(sigma) - (1/(2*sigma^2)) * (log(t) - mu)^2
  
  L_LOG <- sum(status * log_f + (1 - status) * log_S)
  
  return(-L_LOG)
  
}

ajust2 <- optim(par=c(1,1),f=log_vero_surv_lnorm,
                t = tempos,
                status = cens)

###########################

km_c_bx


mu_est <- ajust2$par[1]
sigma_est <- ajust2$par[2]


indices <- which(km_c_bx$sobrevivencia > 0) 
t_km <- km_c_bx$tempo[indices]
s_km <- km_c_bx$sobrevivencia[indices]

# 2. Se o seu km_c_bx não tiver o tempo 0, adicione-o manualmente para o plot:
t_plot <- c(0, t_km)
s_plot <- c(1, s_km)


t_teorico <- seq(0, max(t_plot),length.out=100)
s_lnorm <- pnorm((log(t_teorico) - mu_est) / sigma_est, 
                 lower.tail = FALSE)

plot(t_plot,s_plot, type = "s", lwd = 2,
     xlab = "Tempo (t)", ylab = "S(t)", 
     main = "Ajuste Log-Normal vs. Kaplan-Meier")
lines(t_teorico, s_lnorm, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Kaplan-Meier", "Log-Normal Ajustada"),
       col = c("black", "blue"), lty = c(1, 2), lwd = 2)


# 2. Aplicar a transformação linearizadora
y_coord <- qnorm(1 - s_km) # Inversa da Normal padrão
x_coord <- log(t_km)

# 3. Plotar os dados em torno do modelo
plot(x_coord, y_coord, pch = 16, col = "darkgreen",
     xlab = "log(t)", ylab = "qnorm(1 - S(t))",
     main = "Gráfico de Linearização Log-Normal")

# 4. Adicionar a reta teórica baseada nos seus parâmetros mu e sigma
# y = (1/sigma) * log(t) - (mu/sigma)
abline(a = -mu_est/sigma_est, b = 1/sigma_est, col = "red", lwd = 2)