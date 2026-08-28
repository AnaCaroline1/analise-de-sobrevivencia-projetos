source("turnbull.R")

# =========================================================================
# EXEMPLO 1: SEM GRUPO (curva unica)
# =========================================================================
set.seed(123)
n <- 40
visitas <- seq(2, 24, by = 2)  # pontos de checkup, tipo a tabela da foto

t_verdadeiro <- rweibull(n, shape = 1.5, scale = 12)

fazer_intervalo <- function(t){
  v <- c(0, visitas, Inf)
  k <- findInterval(t, v)
  c(left = v[k], right = v[k + 1])
}

dat <- t(sapply(t_verdadeiro, fazer_intervalo))
data1 <- data.frame(left = dat[, "left"], right = dat[, "right"])

cat("\n===== EXEMPLO 1: sem grupo =====\n")
ajuste1 <- turnbull_boot(data1, B = 200)
print(round(ajuste1, 4))

# comparacao com o pacote survival (mesma NPMLE, so' pra conferir S(t))
ajuste_pkg <- survfit(Surv(left, right, type = "interval2") ~ 1, data = data1)
cat("\nsurvival::survfit (comparacao apenas do S(t)):\n")
print(summary(ajuste_pkg))

# =========================================================================
# EXEMPLO 2: COM GRUPOS (duas curvas, tipo tratamento A vs B)
# =========================================================================
set.seed(456)
n_g <- 40
tA <- rweibull(n_g, shape = 1.5, scale = 10)  # grupo A: eventos mais cedo
tB <- rweibull(n_g, shape = 1.5, scale = 16)  # grupo B: eventos mais tarde

datA <- t(sapply(tA, fazer_intervalo))
datB <- t(sapply(tB, fazer_intervalo))

data2 <- rbind(
  data.frame(left = datA[, "left"], right = datA[, "right"], grupo = "A"),
  data.frame(left = datB[, "left"], right = datB[, "right"], grupo = "B")
)

cat("\n===== EXEMPLO 2: com grupos =====\n")
ajuste2 <- turnbull_by_group(data2, B = 200)
tab2 <- ajuste2[, c("grupo", "time", "n.risk", "n.event", "surv",
                     "std.err", "lower95", "upper95")]
tab2[ , -1] <- round(tab2[ , -1], 4)
print(tab2)

# comparacao com o pacote (estratificado)
ajuste2_pkg <- survfit(Surv(left, right, type = "interval2") ~ grupo, data = data2)
cat("\nsurvival::survfit por grupo (comparacao apenas do S(t)):\n")
print(summary(ajuste2_pkg))

# =========================================================================
# GRAFICO: as duas curvas com IC
# =========================================================================
png("grafico_turnbull_grupos.png", width = 900, height = 600, res = 120)
plot(NULL, xlim = c(0, max(ajuste2$time)), ylim = c(0, 1),
     xlab = "Tempo", ylab = "S(t) estimado (Turnbull)",
     main = "Estimador de Turnbull por grupo")
cores <- c("A" = "steelblue", "B" = "firebrick")
for (g in c("A", "B")){
  sub <- ajuste2[ajuste2$grupo == g, ]
  lines(stepfun(sub$time, c(1, sub$surv)), col = cores[g], lwd = 2, do.points = FALSE)
  lines(stepfun(sub$time, c(1, sub$lower95)), col = cores[g], lty = 2, do.points = FALSE)
  lines(stepfun(sub$time, c(1, sub$upper95)), col = cores[g], lty = 2, do.points = FALSE)
}
legend("topright", legend = c("Grupo A", "Grupo B"), col = cores, lwd = 2, bty = "n")
dev.off()
cat("\nGrafico salvo em grafico_turnbull_grupos.png\n")
