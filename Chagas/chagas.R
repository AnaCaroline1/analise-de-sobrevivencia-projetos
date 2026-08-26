# Pacotes
# ----
library(microdatasus)
library(naniar)
library(lubridate)
library(tidyverse)

# DADOS DO SINAN-SIH-RD - PARTE HOSPITALAR
chagas <- readRDS("sih_rd_23_25.rds") 

colnames(chagas)
table(chagas$CPF_AUT)
df <- chagas |>
  filter(substr(DIAG_PRINC, 1, 3) == "B57") |>
  select(
    # Datas para evolução dos casos:
    DT_INTER, DT_SAIDA, DIAS_PERM, ANO_CMPT,
    # Eventos e Desfecho:
    DIAG_PRINC, COBRANCA, MORTE, CID_MORTE, CID_ASSO,
    # Gravidade, UTI e Procedimentos:
    DIAGSEC1, DIAGSEC2, DIAGSEC3, CAR_INT, PROC_REA,
    UTI_MES_TO, MARCA_UTI, COMPLEX,
    # Covariáveis Demográficas:
    NASC, IDADE, SEXO, RACA_COR, NUM_FILHOS,
    # Geográfica:
    MUNIC_MOV, MUNIC_RES, munResNome, UF_ZI
  ) |>
  mutate(
    across(starts_with("DT_"), ymd),
    ref_data = make_date(year = ANO_CMPT, month = 12, day = 31),
    IDADE = trunc(time_length(interval(NASC, ref_data), unit = "year")),
    # Indicador claro de passagem por UTI (1 = Sim, 0 = Não)
    uso_uti = ifelse(!is.na(UTI_MES_TO) & UTI_MES_TO > 0, 1, 0)
  ) |>
  select(-ref_data)

# dados_brutos <- fetch_datasus(
#     year_start = 2023,
#     year_end = 2025,
#     month_start = 1,
#     month_end = 12,
#     uf = "PA",
#     information_system = "SINAN-CHAGAS"
#   )
# dados_processados <- process_sinan_chagas(dados_brutos)
# write_rds(dados_processados, "chagas_g.rds")
# 
# DADOS DO SINAN-CHAGAS
df_ch_process <- readRDS("chagas_g.rds")
para <- df_ch_process |>
  filter(
    SG_UF == "Pará" | SG_UF == "15",
    SG_UF_NOT == "Pará" | SG_UF_NOT == "15",
    substr(ID_MUNICIP, 1, 2) == "15",
    CLASSI_FIN == "Confirmado" | CLASSI_FIN == 1
  ) |>
  select(
    # Datas para evolução do caso:
    DT_SIN_PRI, DT_NOTIFIC, DT_INVEST, DT_OBITO, DT_ENCERRA,
    # Eventos e Formas de Transmissão:
    CLASSI_FIN, EVOLUCAO, ORAL, ASSINTOMA,
    # Covariáveis Demográficas:
    NU_IDADE_N, CS_SEXO, CS_RACA, CS_ESCOL_N, CS_GESTANT,
    # Covariáveis clínicas (sinais/sintomas de gravidade):
    FEBRE, EDEMA, HEPATOME, ESPLENOM, SINAIS_ICC,
    ARRITMIAS, CHAGOMA, MAECHAGA, MENINGOE, POLIADENO,
    ASTENIA, MICRO_HEMA, OUTRO_SIN,
    # Covariáveis Epidemiológicas/exposição:
    HISTORIA, CON_TRIAT, CON_PROVAV, CRITERIO,
    # Geográfica:
    ID_MUNICIP, SG_UF, SG_UF_NOT, NU_ANO
  ) |>
  mutate(
    unidade_tempo = as.integer(substr(NU_IDADE_N, 1, 1)),
    valor_idade   = as.integer(substr(NU_IDADE_N, 2, 4)),
    idade_anos    = ifelse(unidade_tempo == 4, valor_idade, 0),
    classe_idade  = cut(
      idade_anos,
      breaks = c(-1, 14, 29, 44, 59, 120),
      labels = c("0 a 14 anos", "15 a 29 anos",
                 "30 a 44 anos", "45 a 59 anos", "60 anos ou mais")
    ),
    across(starts_with("DT_"), ymd)
  ) |>
  select(-NU_IDADE_N, -unidade_tempo, -valor_idade)

