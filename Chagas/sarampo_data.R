
# ----
# Dados
# ----
# dados_b <- fetch_datasus(
#   year_start = 2023,
#   year_end = 2025,
#   month_start = 1,
#   month_end = 12,
#   uf = "PA",
#   information_system = "SIH-RD"
# )
# dados_process <- process_sih(dados_b)
# 
# saveRDS(dados_process, "sih_rd_23_25.rds")
# 
sarampo <- readRDS("sih_rd_23_25.rds")

# colnames(sarampo)

# df <- sarampo |>
#   filter(substr(DIAG_PRINC, 1, 3) == "B05") |>
#     select(
#       # Datas para evolução dos casos:
#       DT_INTER,DT_SAIDA,DIAS_PERM,ANO_CMPT,
#       # Eventos:
#       DIAG_PRINC, COBRANCA, MORTE,
#       # Parte clínico e de Gravidade:
#       DIAG_PRINC,DIAG_SECUN, DIAGSEC1,CAR_INT,CAR_INT,PROC_REA,
#       # Covariáveis Demográficas:
#       NASC,IDADE,SEXO,RACA_COR,NUM_FILHOS, DIAS_PERM,
#       # Geográfica:
#       MUNIC_MOV,MUNIC_RES,munResNome
#   )
# saveRDS(df, "sarampo_fill.rds")

# sarampo <- readRDS("sarampo_fill.rds")
# 
# str(sarampo)
# 
# lapply(sarampo, table)
# 
# ----
# Dados Sinan

# library(microdatasus)
# library(naniar)
# library(lubridate)
# library(tidyverse)
# df_sinan <- fetch_datasus(
#   year_start = 2023,
#   year_end = 2025,
#   month_start = 1,
#   month_end = 12,
#   information_system = "SINAN"
# )
# 
# sina_ch <- process_sinan_(df_sinan)
# 
# write_rds(sina_ch, "sinan_chagas")
# 
# dados_sinan <- fetch_datasus(
#   year_start = 2020, 
#   year_end = 2024,
#   information_system = "SINAN"
# )
#----











