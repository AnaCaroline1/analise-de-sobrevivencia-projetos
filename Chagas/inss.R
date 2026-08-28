library(readxl)
library(tidyverse)
library(janitor)
library(lubridate)
# TRATAMENTO BASE DOS DADOS SOBRE AFASTAMENTO POR TRANSTORNOS MENTAIS DO INSS, PEGANDO DO PERÍODO MAIS RECENTE, DE JANEIRO À DEZEMBRO DE 2025.
# 
# A BASE NÃO ESTÁ JUNTA, SÃO 12 PLANILHAS SEPARADAS PARA CADA MÊS. É PRECISO FILTRAR PELO ESTADO E DEPOIS JUNTAR.

inss <- read_excel("GitHub/analise-de-sobrevivencia-projetos/Chagas/DADOS+ABERTOS_CONCEDIDOS_JANEIRO+2024.xlsx", skip = 1, .name_repair = "unique") 

pa <- inss |>
  filter(UF == "Pará",
         Espécie...5 == "Auxílio Doenca Previdenciário",
         # Transtornos mentais de comportamento
         # CID...6 > "F00" & CID...6 < "F99" | CID...6 == "Z73",
         str_detect(CID...6, "^F[0-9]") | str_detect(CID...6, "^Z73")
  ) |>
  mutate(
    `Dt DIB` = ymd(`Dt DIB`),
    `Dt DCB` = case_when(
      `Dt DCB` == "00/00/0000" ~ NA_Date_,
      str_detect(as.character(`Dt DCB`),
                 "^[0-9]+$") ~ excel_numeric_to_date(as.numeric(`Dt DCB`)),
      T ~ as.Date(parse_date_time(`Dt DCB`, orders = c("dmy", "ymd")))),
    status = ifelse(is.na(`Dt DCB`), 0,1),
    dt_lim = dmy("31/01/2024"),
    tempo_dias = as.numeric(coalesce(`Dt DCB`, dt_lim) - `Dt DIB`)
  ) |>
  filter(tempo_dias > 0) |>
  select(-`País de Acordo Internacional`)
