log_verossimilhança_weibull <- function(theta,t){
  gama <- theta[1]
  alpha <- theta[2] #gamma #alpha
  
  # Evitando parâmetros não-positivos durante a otimização
  if(gama <= 0 || alpha <=0) return(Inf)
  
  # Log da densidade weibull
  
  log_link_w <- dweibull(t,shape=gama,scale=alpha,log = T)
  
  return(-sum(log_link_w))
}

# Estimando parâmetros a partir dos meus dados como a função fitdistr do pacote MASS

tempos # Do banco de dados sobre câncer de bexiga

ajust_ini1 <- fitdistr(tempos, "weibull")

par_weibull1 <- ajust_ini$estimate


log_verossimilhança_weibull(par_weibull,tempos)

# Aplicando a optimização

fit1 <- optim(par=c(1,1),fn=log_verossimilhança_weibull, t = tempos)

gama_est1 <- fit1$par[1]
alph_est1 <- fit2$par[2]

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
?lines
