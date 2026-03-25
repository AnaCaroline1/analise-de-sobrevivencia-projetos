# Nova função considerando censura
log_vero_surv_exponencial <- function(theta, t, status) {
  alpha <- theta[1]
  
  if(alpha <= 0) return(Inf)
  
  logv <- sum(-status * log(alpha) -(t/alpha))
  
  # Log-verossimilhança total (negativa para o optim)
  return(-logv)
}

# Estimando parâmetros a partir dos meus dados como a função fitdistr do pacote MASS

tempos # Do banco de dados sobre câncer de bexiga
cens
ajust_ini1 <- fitdistr(tempos, "exponential")

par_exp1 <- ajust_ini1$estimate


log_vero_surv_exponencial(par_exp1, tempos, cens)


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
                method="Brent",
                lower = 0.001,
                upper = 1000,
                hessian = T)





n <- length(tempos) # tamanho da amostra
x <- rexp(n,ajus11$par)

x_seq <- seq(min(tempos), max(tempos), length.out = 10)

#length.out: gera uma sequência de 100 valores x igualmente espaçados entre o valor
##mínimo e o valor máximo dos dados.
##
densidade_exp <- dexp(x_seq,ajus11$par) 

plot(x, ylab = "Frequência", xlab = "Valores dos dados", main = "Histograma com ajuste: Exponencial p/ câncer de bexiga",
     freq = FALSE)
lines(x_seq, densidade_exp, col = "red", lwd = 2)
hist(x, ylab = "Frequência", xlab = "Valores dos dados", main = "Histograma com ajuste: Exponencial p/ câncer de bexiga",
     freq = FALSE)
lines(x_seq, densidade_exp, col = "red", lwd = 2)

#############################################
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

km_c_bx <- ekm(tempos,cens)
km_c_bx

t_teorico <- seq(0, max(tempos),length.out=20)
s_teorico <- exp(-t_teorico/ajus11$par)

plot(km_c_bx$tempo, km_c_bx$sobrevivencia,
     type = "s",
     lwd=2,
     main = "Ajuste Manual: Dados vs. Modelo Exponencial",
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


#############################################

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

###########################################

# Criando intervalo de confiança para os resultados obtidos

# Pegando a hessiana dos ajuste

f_theta <- ajus11$hessian


# Variância e EP de cada distribuição
# Calculando a Var(theta) 

var_theta <- f_theta

ep_theta <- 

## Calculando os intervalos de confiança

limite_inf <- round(exp(-(ajus11$par - 1.96 * as.vector(ep_theta)))*t_teorico,4)
limite_sup <- round(exp(-(ajus11$par + 1.96 * as.vector(ep_theta)))*t_teorico,4)

## Plotando em gráfico com intervalo de confiança
 

# Criando o gráfico

c(t_teorico, rev(t_teorico)) == c(limite_sup, rev(limite_inf))

plot(km_c_bx$tempo, km_c_bx$sobrevivencia,
     type = "n",
     lwd=2,
     main = "Ajuste Manual: Dados vs. Modelo Exponencial",
     ylim = c(min(limite_inf),max(limite_sup)),
     xlab = "Tempo", ylab = "Sobrevivência"
     )
polygon(c(t_teorico, rev(t_teorico)),
        c(limite_sup, rev(limite_inf)),
        col = rgb(1,0,0, alpha = 0.15),
        border = NA
        )
lines(km_c_bx$tempo, km_c_bx$sobrevivencia, 
      type = "s", lwd = 1)
lines(t_teorico, s_teorico, col = "red", lwd = 2, lty = 2)
legend("topright", 
       legend = c("Kaplan-Meier (Dados)", "Exponencial (Modelo)", "IC 95% (modelo)"),
       col = c("black", "red", rgb(1,0,0,alpha=0.3)), 
       lty = c(1, 2, 1), 
       lwd = c(2,2,8)
       )

#ggplot(tempos, aes(x = grupo, y = media)) +
#  geom_bar(stat = "identity", fill = "skyblue") +
#  geom_errorbar(aes(ymin = limite_inf, ymax = limite_sup), width = 0.2) +
#  theme_minimal() 



###########################################


