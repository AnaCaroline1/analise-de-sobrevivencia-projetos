library(readxl)
library(dplyr)
library(tidyverse)
library(janitor)
library(janitor) # Formatação de datas
library(arrow) # Formato mais rápido de leitura
library(purrr)
library(stringr)
# Tratamento das bases sobre afastamento por transtornos mentais do INSS, ano de 2025, o mais recente e completo com os 12 meses;
# ----
# O objetivo aqui é juntar as bases já que é relativo a um mês por base, então de forma sistemática o que quero aplicar:
# 1. Carregar os bancos;
# 2. Filtrar pelo estado do Pará;
# 3. Em espécie (Tipo de benefício consedido), apenas os "Auxílio Doenca Previdenciário", que são de nosso interesse;
# 4. Aplicar tratamento para as datas (pensando em aplicar por último, depois de filtrar nos passos 1-3);
# 5. Por fim filtrar os CID's, str_detect(CID...6, "^F[0-9]") | str_detect(CID...6, "^Z73") -> Esses são os cid's que são por afastamento por Transtornos Mentais (ansiedade e seus quadros, depressão e seus quados etc.), como minha referência é da base de 2024, pode ser que tenha (espero), mais observações para burnout e correlacionados, isto será verificado ao decorrer dos filtros.
# 
# Pensando bem, poderia deixar os cid's por último, apenas aplicando as etapas de 1-3, depois filtrando os cids, conforme a importância para o estudo.
# 
# Dic. Nomes
# ----

# ----
# jan <- read_excel("INSS_TM/CONCEDIDOS_DADOS_ABERTOS_JANEIRO+2025.xlsx", skip = 1, .name_repair = "unique")
# 
# 
# jan_1 <- jan |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#          )
# write_parquet(jan_1, "jan_fit.parquet")
         
# ----
# fev <- read_excel("INSS_TM/CONCEDIDOS_DADOS+ABERTOS_FEVEREIRO+2025.xlsx", skip = 1, .name_repair = "unique")
# 
# fev_1 <-  fev |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#   )
# 
# write_parquet(fev_1, "fev_fit.parquet")

# ----
# mar <- read_excel("INSS_TM/BEN_CONCEDIDOS_032025.xlsx", skip = 1, .name_repair = "unique")

# mar_1 <- mar |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#   )
# write_parquet(mar_1, "mar_fit.parquet")
# ----

# abr <- read_excel("INSS_TM/BEN_CONCEDIDOS_042025.xlsx", skip = 1, .name_repair = "unique")
# 
# abr_1 <- abr |> 
#   filter(UF == "Pará",
#        Espécie...5 == "Auxílio Doenca Previdenciário",
#        # Transtornos mentais e comportamentais
#        str_detect(CID...6, "^F[0-9]") |
#          # Síndrome de Burnout
#          str_detect(CID...6, "^Z73")
#   )
# write_parquet(abr_1, "abr_fit.parquet")
# ----
# mai <- read_excel("INSS_TM/CONCEDIDOS+DADOS+ABERTOS+MAIO+2025.xlsx", skip = 1, .name_repair = "unique")
# mai_1 <- mai |>
#   filter(UF == "Pará",
#        Espécie...5 == "Auxílio Doenca Previdenciário",
#        # Transtornos mentais e comportamentais
#        str_detect(CID...6, "^F[0-9]") |
#          # Síndrome de Burnout
#          str_detect(CID...6, "^Z73")
#   )
# write_parquet(mai_1, "mai_fit.parquet")
# ----

# jun <- read_excel("INSS_TM/CONCEDIDOS+DADOS+ABERTOS+JUNHO+2025.xlsx", skip = 1, .name_repair = "unique")
# 
# jun_1 <- jun |>
#   filter(UF == "Pará",
#        Espécie...5 == "Auxílio Doenca Previdenciário",
#        # Transtornos mentais e comportamentais
#        str_detect(CID...6, "^F[0-9]") |
#          # Síndrome de Burnout
#          str_detect(CID...6, "^Z73")
#   )
# write_parquet(jun_1, "jun_fit.parquet")
# ----
jul <- read_excel("INSS_TM/CONCEDIDOS+DADOS+ABERTOS_JULHO+2025.xlsx")
jul_1 <- jul |>
  filter(UF == "Pará",
         Espécie...5 == "Auxílio Doenca Previdenciário",
         # Transtornos mentais e comportamentais
         str_detect(CID...6, "^F[0-9]") |
           # Síndrome de Burnout
           str_detect(CID...6, "^Z73")
  )
write_parquet(jul_1, "jul_fit.parquet")
# ----

# ago <- read_excel("INSS_TM/CONCEDIDOS_DADOS+ABERTOS_AGOSTO+DE+2025.xlsx")
# ago_1 <- ago |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#   )
# write_parquet(ago_1, "ago_fit.parquet")
# ----
# set <- read_excel("INSS_TM/CONCEDIDOS+DADOS+ABERTOS_SETEMBRO+2025.xlsx", skip = 1, .name_repair = "unique")
# 
# set_1 <- set |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#   )
# write_parquet(set_1, "set_fit.parquet")

# ----
# Outubro precisa ser tratada uma vez que em cid6 não foi preenchido os cid, mas se tem indicação de cid em cid 7, vai dar um pouco mais de trabalho.
# ----
# out <- read_excel("INSS_TM/CONCEDIDOS+DADOS+ABERTOS_OUTUBRO+2025.xlsx", skip = 1, .name_repair = "unique")
# 
# out_1 <- out |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário"
#          # # Transtornos mentais e comportamentais
#          # str_detect(CID...7, "^F"),
#          # # Síndrome de Burnout
#          # str_detect(CID...7, "^Z")
#   )
# write_parquet(out_1, "out_fit.parquet")
# ----
# nov <- read_excel("INSS_TM/DADOS+ABERTOS_CONCEDIDOS_NOVEMBRO+2025.xlsx", skip = 1, .name_repair = "unique")
# 
# nov_1 <- nov |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#   )
# write_parquet(nov_1, "nov_fit.parquet")

# dez <- read_excel("INSS_TM/BEN_CONCEDIDOS_122025.xlsx", skip = 1, .name_repair = "unique")
# 
# dez_1 <- dez |>
#   filter(UF == "Pará",
#          Espécie...5 == "Auxílio Doenca Previdenciário",
#          # Transtornos mentais e comportamentais
#          str_detect(CID...6, "^F[0-9]") |
#            # Síndrome de Burnout
#            str_detect(CID...6, "^Z73")
#   )
# write_parquet(dez_1, "dez_fit.parquet")


# ----
# Juntanto as bases
# ----
# arq <- list.files(path = "fit_inss/", pattern = "\\.parquet$", full.names = T)
# 
# limpar_nomes <- function(df) {
#   novos_nomes <- names(df) %>%
#     str_remove("\\.\\.\\.\\d+$")
#   names(df) <- make.unique(novos_nomes, sep = "_")
#   df
# }
# 
# lista_bases <- map(
#   arq,
#   ~ read_parquet(.x) %>%
#     limpar_nomes() %>%
#     mutate(across(everything(), as.character))   # força tudo pra texto
# )
# 
# inss_pa <- bind_rows(lista_bases)

#write_parquet(inss_pa, "inss_pa.parquet")


# inss_pa <- read_parquet("inss_pa.parquet")
# 
# 
# inss_out <- read_parquet("out_fit.parquet")
# 
# fit_t <- inss_pa |>
#   mutate(
#       CID = if_else(
#         `Competência concessão` == "202510" & (is.na(CID) | CID == ""),
#         str_extract(str_trim(CID_1), "^[A-Z][0-9]{2}(\\.[0-9]{1,2})?"),
#         CID
#       )
#   ) |>
#     filter(
#            # Transtornos mentais e comportamentais
#            str_detect(CID, "^F[0-9]") |
#              # Síndrome de Burnout
#              str_detect(CID, "^Z73")
#     )
  
#write_parquet(fit_t, "inss_pa.parquet")
