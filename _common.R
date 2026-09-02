library(ggplot2)
library(tidyverse)
library(janitor)
library(knitr)
# Precisa do pacote arrow
#install.packages("arrow")
library(arrow)
library(readxl)
library(MASS)
# Carregar pacotes
library(survival)
library(survminer) # Para gráficos (opcional)

tempos<-c(3,5,6,7,8,9,10,10,12,15,15,18,19,20,22,25,28,30,40,45)
cens<-c(1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0)

df_cc_bx <- data.frame(tempos, cens)

set.seed(42)

df_cc_bx$grupos <- sample(c("A","B"),size = nrow(df_cc_bx),replace = T)
df_cc_bx$idade <- round(runif(nrow(df_cc_bx),min = 20, max = 60))

df_cc_bx
#write_parquet(df_final, "df_srag_sg.parquet")
# 
# dengue <- read_parquet("df_dengue_2025.parquet")
