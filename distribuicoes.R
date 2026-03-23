# Funções de log-verossimilhança aplicadas

log_vero_surv_exponencial <- function(theta, t, status) {
  alpha <- theta[1]
  
  if(alpha <= 0) return(Inf)
  
  logv <- sum(-status * log(alpha) -(t/alpha))
  
  # Log-verossimilhança total (negativa para o optim)
  return(-logv)
}

log_vero_surv_lnorm <- function(theta,t,status){
  mu    <- theta[1]
  sigma <- theta[2]
  log_S <- plnorm(t, meanlog = mu, sdlog = sigma, lower.tail = FALSE, log.p = TRUE)
  
  log_f <- -0.5*log(2*pi) - log(t) - log(sigma) - (1/(2*sigma^2)) * (log(t) - mu)^2
  
  L_LOG <- sum(status * log_f + (1 - status) * log_S)
  
  return(-L_LOG)

}

log_vero_surv_weibull <- function(theta, t, status){
  
  alpha <- theta[1] 
  gamma <- theta[2]
  r <- status
  
  Log_L <- sum(r*(log(gamma)-gamma*log(alpha)+(gamma-1)*log(t))-(t/alpha)^gamma)
  
  if(!is.finite(Log_L)) return(1e10) #garante o tempo positivo apenas
  
  return(-Log_L)
}
S_loglogistica <- function(theta, t, delta) {
  
  alpha <- exp(theta[1])
  gamma <- exp(theta[2])
  
  log_f <- log(gamma) - (gamma * log(alpha)) + (gamma - 1) * log(t) - 2 * log(1 + (t/alpha)^gamma)
  
  log_S <- -log(1 + (t/alpha)^gamma)
  
  LOG_L <- sum(delta * log_f + (1 - delta) * log_S)
  
  return(-LOG_L)
}
###########################################################

# Estimando parâmetros a partir dos meus dados como a função fitdistr do pacote MASS

tempos # Do banco de dados sobre câncer de bexiga
cens
ajust_ini1 <- fitdistr(tempos, "exponential")

par_exp1 <- ajust_ini1$estimate

###########################################################

# Aplicando a optimização

ajus1 <- optimize(f=log_vero_surv_exponencial, 
                  interval = c(0.001, 1000),
                  t = tempos,
                  status = cens)

ajus1$minimum

ajus11 <- optim(par = 1,
                fn=log_vero_surv_exponencial, 
                t = tempos,
                status = cens, 
                method="Brent",#para um parâmetro
                lower = 0.001,
                upper = 1000)
alpha_est <- ajus11$par

ajust2 <- optim(par=c(1,1),f=log_vero_surv_lnorm,
                   t = tempos,
                   status = cens)
ajust3 <- optim(par=c(1,1),f=log_vero_surv_weibull,
                t = tempos,
                status = cens)

# Chute inicial (log da média dos tempos e log de 1)
inicio <- c(log(mean(tempos)), 0)

fit_ll <- optim(par = inicio, fn = S_loglogistica, t = tempos, 
                delta = cens, method = "BFGS")
fit_ll
############################################################

## Estimador KM
## Passos
# ordenar os tempos
# calcular a proporção de sobreviventes para
# o tempo de evento (delta=1)

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
##########################################################

km_c_bx <- ekm(tempos,cens)
km_c_bx
#########################################################
# EXPONENCIALLLLLL
t_teorico <- seq(0, max(tempos),length.out=100)
s_teorico <- exp(-t_teorico/ajus11$par)

plot(km_c_bx$tempo, km_c_bx$sobrevivencia,
     type = "s",
     lwd=2,
     main = "Ajuste Manual: Dados(pacientes com câncer de bexiga) vs. Modelo Exponencial",
     ylim = c(0,1))
lines(t_teorico,s_teorico, col="red",lwd=2,lty=2)
legend("topright", 
       legend = c("Kaplan-Meier (Dados)", "Exponencial (Modelo)"),
       col = c("black", "red"), 
       lty = c(1, 2), 
       lwd = 2)

# Para entender se a exponencial é um boa escolha, se utiliza
# uma visualização gráfica por meio do Risco Acumulado, sendo na 
# exponencial -log(S(t)), com linha reta começando em zero.

# Calcular o log negativo da sobrevivência de Kaplan-Meier
risco_acumulado_km <- -log(km_c_bx$sobrevivencia)

# Plotar
plot(km_c_bx$tempo, risco_acumulado_km, pch = 5, col = "blue",
     xlab = "Tempo (t)", ylab = "-log(S(t))",
     main = "Verificação de Risco Constante")

# Adicionar a linha esperada pelo seu modelo (1/alpha * t)
abline(a = 0, b = 1/ajus11$par, col = "red", lwd = 2)

# Se os pontos azuis seguirem a reta, a 
#distribuição é adequada para os dados.
# Se os pontos fazem uma curva, seja para cima ou baixo,
# indica que o risco não é constante, logo se testa outra 
#distribuição.
#######################################################
# Log normal

mu_est <- ajust2$par[1]
sigma_est <- ajust2$par[2]

s_lnorm <- pnorm((log(t_teorico) - mu_est) / sigma_est, 
                 lower.tail = FALSE)

plot(km_c_bx$tempo, km_c_bx$sobrevivencia, type = "s", lwd = 2,
     xlab = "Tempo (t)", ylab = "S(t)", 
     main = "Ajuste Log-Normal vs. Kaplan-Meier")
lines(t_teorico, s_lnorm, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = c("Kaplan-Meier", "Log-Normal Ajustada"),
       col = c("black", "blue"), lty = c(1, 2), lwd = 2)

# 1. Pegar apenas os tempos onde houve evento (evita erros no qnorm)
indices <- which(km_c_bx$sobrevivencia > 0 & km_c_bx$sobrevivencia < 1)
t_km <- km_c_bx$tempo[indices]
s_km <- km_c_bx$sobrevivencia[indices]

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
#############################################
# Weibull


alpha_est <- ajust3$par[1]
gamma_est <- ajust3$par[2]
t_seq <- seq(0, max(tempos), length.out = 100)
s_weibull <- exp(-(t_seq / alpha_est)^gamma_est)

# 3. Plotar
plot(km_c_bx$tempo, km_c_bx$sobrevivencia, type = "s", lwd = 2,
     xlab = "Tempo (t)", ylab = "S(t)", 
     main = "Ajuste Weibull vs. Kaplan-Meier")
lines(t_seq, s_weibull, col = "darkred", lwd = 2, lty = 2)
legend("topright", legend = c("Kaplan-Meier", "Weibull Ajustada"),
       col = c("black", "darkred"), lty = c(1, 2), lwd = 2)


# 1. Preparar os dados (remover S(t)=0 ou 1 para evitar log de erro)
indices <- which(km_c_bx$sobrevivencia > 0 & km_c_bx$sobrevivencia < 1)
t_km <- km_c_bx$tempo[indices]
s_km <- km_c_bx$sobrevivencia[indices]

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
###############################################
#Log_logistica
#
t_seq <- seq(0.0001, max(tempos, na.rm = TRUE), length.out = 500)

# 2. Calcular a sobrevivência teórica (Log-Logística)
# S(t) = 1 / (1 + (t/alpha)^gamma)
s_loglogis <- 1 / (1 + (t_seq / alpha_est)^gamma_est)

# 3. Plotagem com tratamento de erro
# Usei km_c_bx$tempo e km_c_bx$sobrevivencia como você definiu
plot(km_c_bx$tempo, km_c_bx$sobrevivencia, 
     type = "s", lwd = 2,
     xlab = "Tempo (t)", ylab = "S(t)", 
     main = "Ajuste Log-Logístico vs. Kaplan-Meier",
     ylim = c(0, 1)) # Garante que o eixo Y vá de 0 a 1

# Adiciona a linha do modelo
lines(t_seq, s_loglogis, col = "darkgreen", lwd = 2, lty = 2)

# Adiciona a legenda
legend("topright", 
       legend = c("Kaplan-Meier", "Log-Logística Ajustada"),
       col = c("black", "darkgreen"), 
       lty = c(1, 2), 
       lwd = 2,
       bty = "n") # bty="n" remove a caixa
###############################################
par_est_1 <- fit1$par
shape_final <- par_est_1[1]
scale_final <- par_est_1[2]
s_weibull <- function(t) {
  p_shape <- par_est_1[1]
  p_scale <- par_est_1[2]
  return(exp(-(t / p_scale)^p_shape))
}

# Plotando
tempos <- seq(0, max(tempos), length.out = 100)
plot(tempos, s_weibull(tempos), type = "l", col = "blue", lwd = 2,
     main = "Curva de Sobrevivência Estimada (Weibull)",
     xlab = "Tempo", ylab = "S(t)")

S_t <- function(tempo_alvo) {
  # pweibull dá a área à esquerda (F(t)), 
  # então 1 - pweibull dá a área à direita (Sobrevivência)
  pweibull(tempo_alvo, shape = shape_final, scale = scale_final, lower.tail = FALSE)
}

S_t(30)

##############################################

# Verificando como os dados do livro: câncer de bexiga

# Ajuste conforme a distribuição

ajust1_cc_bx <- survreg(Surv(tempos,cens)~1,dist = "exponential")
ajust2_cc_bx <- survreg(Surv(tempos,cens)~1,dist = "weibull")
ajust3_cc_bx <- survreg(Surv(tempos,cens)~1,dist = "lognorm")
ajust4_cc_bx <- survreg(Surv(tempos,cens)~1,dist = "loglogistic")
# Ekm  dos dados

ekm_cc_bx <- survfit(Surv(tempos,cens)~1)
time <- ekm_cc_bx$time
st <- ekm_cc_bx$surv

# Cálculo dos parâmetros estimados das distribuições

alpha<- round(exp(ajust1_cc_bx$coefficients[1]),2 )
alpha2<-  round(exp(ajust2_cc_bx$coefficients[1]),2 )
gama2<-1/ajust2_cc_bx$scale

mi <-  round(ajust3_cc_bx$coefficients[1],2)
alpha3 <- round(ajust3_cc_bx$scale,2)

ste_cc_bx<-exp(-time/alpha)

stw_cc_bx<- exp(-(time/alpha2)^gama2)

stln_cc_bx<- pnorm((-log(time)+ mi)/alpha3)

#par(mfrow=c(1,3))
plot(st, ste_cc_bx, pch=16, ylim=c(0,1), xlim=c(0,1), 
     xlab = "S(t): Kaplan-Meier", ylab="S(t): exponencial",
     cex.lab=1.5,cex.axis=1.3)
abline(a=0,b=1)

plot(st, stw_cc_bx, pch=16, ylim=c(0,1), xlim=c(0,1), 
     xlab = "S(t): Kaplan-Meier", ylab="S(t): weibull",
     cex.lab=1.5,cex.axis=1.3)
abline(a=0,b=1)

plot(st, stln_cc_bx, pch=16, ylim=c(0,1), xlim=c(0,1), 
     xlab = "S(t): Kaplan-Meier", ylab="S(t): log-normal",
     cex.lab=1.5,cex.axis=1.3)
abline(a=0,b=1)

par(mfrow=c(1,3),mar=c(5,5,2,2)); 
invst<-qnorm(st) 
plot(time,-log(st),pch=16,xlab="Tempo",ylab="-log(S(t))") 
plot(log(time),log(-log(st)),pch=16,xlab="log(Tempo)",ylab="log(-log(S(t)))") 
plot(log(time),invst,pch=16,xlab="log(Tempo)",ylab=expression(Phi^-1*(S(t))))

par(mfrow=c(1,2)) 
plot(ekm_cc_bx,conf.int=F,xlab="Tempo",ylab="S(t)"); lines(c(0,time),c(1,stw_cc_bx),lty=2) 
legend(25,0.8,lty=c(1,2),c("Kaplan-Meier","Weibull"),bty="n",cex=0.9) 
plot(ekm_cc_bx,conf.int=F,xlab="Tempo",ylab="S(t)"); lines(c(0,time),c(1,stln_cc_bx),lty=2) 
legend(25,0.8,lty=c(1,2),c("Kaplan-Meier","Log-normal"),bty="n",cex=0.9)




