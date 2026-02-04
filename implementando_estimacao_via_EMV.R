log_verossimilhança_weibull <- function(theta,t){
  gama <- theta[1]
  alpha <- theta[2] #gamma #alpha
  
  # Evitando parâmetros não-positivos durante a otimização
  if(gama <= 0 || alpha <=0) return(Inf)
  
  # Log da densidade weibull
  
  log_link_w <- dweibull(t,shape=gama,scale=alpha,log = T)
  
  return(-sum(log_link_w))
}

# Nova função considerando censura
log_verossimilhanca_censura <- function(theta, t, status) {
  gama <- theta[1]  # shape
  alpha <- theta[2] # scale
  
  if(gama <= 0 || alpha <= 0) return(Inf)
  
  # Evento de falha (Densidade)
  # dweibull(..., log = T)
  f_t <- dweibull(t[status == 1], shape = gama, scale = alpha, log = TRUE)
  
  # Evento de censura (Sobrevivência)
  # pweibull(..., lower.tail = F, log = T) é o log(S(t))
  
  S_t <- pweibull(t[status == 0], shape = gama, scale = alpha, lower.tail = FALSE, log.p = TRUE)
  
  # Log-verossimilhança total (negativa para o optim)
  return(-(sum(f_t) + sum(S_t)))
}

# Estimando parâmetros a partir dos meus dados como a função fitdistr do pacote MASS

tempos # Do banco de dados sobre câncer de bexiga
cens
ajust_ini1 <- fitdistr(tempos, "weibull")

par_weibull1 <- ajust_ini1$estimate


log_verossimilhança_weibull(par_weibull1,tempos)

# Aplicando a optimização

fit1 <- optim(par=c(1,1),fn=log_verossimilhanca_censura, t = tempos,status = cens)

gama_est1 <- fit1$par[1]
alph_est1 <- fit1$par[2]

n <- length(tempos) # tamanho da amostra
x <- rweibull(n, shape=gama_est1, scale = alph_est1)

x_seq <- seq(min(tempos), max(tempos), length.out = 10)
#length.out: gera uma sequência de 100 valores x igualmente espaçados entre o valor
##mínimo e o valor máximo dos dados.
densidade_weibull1 <- dweibull(x_seq, shape = gama_est1, scale = alph_est1) 

hist(x, ylab = "Frequência", xlab = "Valores dos dados", main = "Histograma com ajuste: Weibull p/ bexiga",
     freq = FALSE)
lines(x_seq, densidade_weibull1, col = "red", lwd = 2)

# Estimando parâmetros com os bancos de dados.

tempo_sg # Do banco de dados sobre Sg e Srag

ajust2 <- fitdistr(tempo_sg, "weibull")
par_weibull2 <- ajust2$estimate

log_verossimilhança_weibull(par_weibull2,tempo_sg)

fit2 <- optim(par=c(1,1),fn=log_verossimilhança_weibull, t = tempo_sg)


gama_est <- fit2$par[1]
alph_est <- fit2$par[2]

n <- length(tempo_sg) # tamanho da amostra
x <- rweibull(n, shape=gama_est, scale = alph_est)

x_seq <- seq(min(tempo_sg), max(tempo_sg), length.out = 5)
#length.out: gera uma sequência de 100 valores x igualmente espaçados entre o valor
##mínimo e o valor máximo dos dados.
densidade_weibull2 <- dweibull(x_seq, shape = gama_est, scale = alph_est) 

hist(x, ylab = "Frequência", xlab = "Valores dos dados", main = "Histograma com ajuste: Weibull p/ srg",
     freq = FALSE)
lines(x_seq, densidade_weibull2, col = "red", lwd = 2)

## Estimador para os dados sobre câncer de pele

df_cc_pl

tempo_3 <- df_cc_pl$survtime
cens_3 <- df_cc_pl$status
ajus3 <- fitdistr(tempo_3,"weibull")
par_weibull3 <- ajus3$estimate
log_verossimilhanca_censura(par_weibull3,tempo_3,cens_3)

fit3 <- optim(par = c(1,1),fn=log_verossimilhanca_censura,t=tempo_3,status = cens_3)



gama_est3 <- fit3$par[1]
alph_est3 <- fit3$par[2]

n <- length(tempo_3) # tamanho da amostra
x <- rweibull(n, shape=gama_est3, scale = alph_est3)

x_seq <- seq(min(tempo_3), max(tempo_3), length.out = 100)
#length.out: gera uma sequência de 100 valores x igualmente espaçados entre o valor
##mínimo e o valor máximo dos dados.
densidade_weibull3 <- dweibull(x_seq, shape = gama_est3, scale = alph_est3) 


par(mar=c(5,4,4,2))
hist(x, ylab = "Frequência", xlab = "Valores dos dados", main = "Histograma com ajuste: Weibull p/ pele",
     freq = FALSE)
lines(x_seq, densidade_weibull3, col = "red", lwd = 2)


#############################################

par_est_1 <- fit1$par

par_est_3 <- fit3$par
shape_final <- par_est_1[1]
scale_final <- par_est_1[2]

s_weibull <- function(par_f,t) {
  p_shape <- par_f[1]
  p_scale <- par_f[2]
  return(exp(-(t / p_scale)^p_shape))
}
s_weibull(par_est_3,tempo_3)
# Plotando
tempos <- seq(0, max(tempos), length.out = 100)
plot(tempos, s_weibull(tempos), type = "l", col = "blue", lwd = 2,
     main = "Curva de Sobrevivência Estimada (Weibull)",
     xlab = "Tempo", ylab = "S(t)")

tempo3 <- seq(0, max(tempo_3), length.out = 100)

plot(tempo3, s_weibull(par_est_3,tempo_3), type = "l", col = "blue", lwd = 2,
     main = "Curva de Sobrevivência Estimada (Weibull)",
     xlab = "Tempo", ylab = "S(t)")

S_t <- function(tempo_alvo) {
  # pweibull dá a área à esquerda (F(t)), 
  # então 1 - pweibull dá a área à direita (Sobrevivência)
  pweibull(tempo_alvo, shape = shape_final, scale = scale_final, lower.tail = FALSE)
}

S_t(30)
