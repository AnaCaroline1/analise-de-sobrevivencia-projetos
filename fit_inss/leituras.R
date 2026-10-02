library(arrow)

library(janitor)
library(digest)
library(dplyr)

library(tidyverse)
library(plotly)
library(scales)

source("inss_sobrevivencia.R")

cov_data <- function(x) {
  x <- trimws(x)
  case_when(
    grepl("^[0-9]{4,6}$", x) & as.numeric(x) < 60000 ~
      janitor::excel_numeric_to_date(as.numeric(x)),
    # Qualquer outro formato de data textual: tenta múltiplas ordens
    TRUE ~ as.Date(lubridate::parse_date_time(x, orders = c(
      "ymd", "dmy", "mdy", "Ymd"
    )))
  )
}
inss <- read_parquet("inss_pa.parquet") |>
  clean_names() |>
  mutate(
    across(starts_with("dt_"), ~ cov_data(.x)),
    qt_sm_rmi = as.numeric(qt_sm_rmi),
    qt_anos_contribuicao = as.numeric(qt_anos_contribuicao)
  )

# --- 1. Prep: variáveis derivadas de tempo ---
inss <- inss %>%
  mutate(
    idade_concessao = as.numeric(difftime(dt_dib, dt_nascimento, units = "days")) / 365.25,
    duracao_beneficio_dias = as.numeric(difftime(
      coalesce(dt_dcb, Sys.Date()), dt_dib, units = "days"
    )),
    status_beneficio = if_else(is.na(dt_dcb), "Ativo (não cessado)", "Cessado"),
    ano_mes_dib = floor_date(dt_dib, "month")
  )

# --- 2. CID mais frequentes ---
top_cid <- inss %>%
  count(cid_1, sort = TRUE) %>%
  slice_max(n, n = 15) %>%
  mutate(cid_1 = fct_reorder(cid_1, n))

p_cid <- ggplot(top_cid, aes(
  x = n,
  y = cid_1,
  text = paste0(cid_1, ": ", n)
)) +
  geom_col(fill = "#2c7fb8") +
  labs(title = "CIDs mais frequentes", x = "Nº de benefícios", y = NULL) +
  theme_minimal()


# --- 3. Distribuição por sexo ---
p_sexo <- inss %>%
  count(sexo) %>%
  plot_ly(
    labels = ~ sexo,
    values = ~ n,
    type = "pie",
    textinfo = "label+percent",
    marker = list(colors = c("#2c7fb8", "#f03b20"))
  ) %>%
  layout(title = "Distribuição por sexo")
p_sexo

# --- 4. Idade na concessão do benefício ---
p_idade <- ggplot(inss, aes(x = idade_concessao, fill = sexo)) +
  geom_histogram(binwidth = 5,
                 alpha = 0.7,
                 position = "identity") +
  labs(title = "Idade na concessão do benefício", x = "Idade (anos)", y = "Frequência") +
  theme_minimal()



# --- 6. Duração por CID (top 10) — boxplot ---
top10_cid_codes <- inss %>% count(cid, sort = TRUE) %>% slice_max(n, n = 10) %>% pull(cid)

p_duracao_cid <- inss %>%
  filter(cid %in% top10_cid_codes) %>%
  ggplot(aes(
    x = fct_reorder(cid, duracao_beneficio_dias, .fun = median),
    y = duracao_beneficio_dias / 30
  )) +
  geom_boxplot(fill = "#41b6c4") +
  coord_flip() +
  labs(title = "Duração do benefício por CID (top 10)", x = "CID", y = "Meses") +
  theme_minimal()



# --- 7. Salário (RMI em salários mínimos) ---
p_salario <- ggplot(inss, aes(x = qt_sm_rmi)) +
  geom_histogram(bins = 30, fill = "#238443") +
  labs(title = "Renda Mensal Inicial (em salários mínimos)", x = "Qtd. salários mínimos", y = "Frequência") +
  theme_minimal()



# --- 8. Ramo de atividade mais frequente ---
top_ramo <- inss %>%
  count(ramo_atividade, sort = TRUE) %>%
  slice_max(n, n = 10) %>%
  mutate(ramo_atividade = fct_reorder(ramo_atividade, n))

p_ramo <- ggplot(top_ramo, aes(x = n, y = ramo_atividade)) +
  geom_col(fill = "#fd8d3c") +
  labs(title = "Ramo de atividade (top 10)", x = "Nº de benefícios", y = NULL) +
  theme_minimal()


# --- 9. Geografia: top municípios ---
top_mun <- inss %>%
  count(mun_resid, sort = TRUE) %>%
  slice_max(n, n = 15) %>%
  mutate(mun_resid = fct_reorder(mun_resid, n))

p_mun <- ggplot(top_mun, aes(x = n, y = mun_resid)) +
  geom_col(fill = "#756bb1") +
  labs(title = "Municípios com mais concessões (top 15)", x = "Nº de benefícios", y = NULL) +
  theme_minimal()

# --- 10. Série temporal: concessões ao longo do tempo ---
p_serie <- inss %>%
  count(ano_mes_dib) %>%
  ggplot(aes(x = ano_mes_dib, y = n)) +
  geom_line(color = "#2c7fb8") +
  geom_point(color = "#2c7fb8", size = 1) +
  labs(title = "Concessões de benefício ao longo do tempo", x = "Mês/Ano (DIB)", y = "Nº de concessões") +
  scale_x_date(date_labels = "%b/%Y", date_breaks = "3 months") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ----
# Identificando o mesmo indivíduo, se existir, na base e identificar se tecnicamente são duplicados os dados ou se o mesmo indivíduo recebeu benefícios diferentes
# ----
# Pseudo-ID
# Como funciona: combina variáveis identificadoras que possuem alta discriminação
# ----
# No banco dá para combinar a data de nascimento, município/cep de residência, gênero, escolariadade (essa última não tem muitas células preenchidas, então não é tão utilizável para a proposta);
# ----

inss_ps <- inss |>
  rowwise() |>
  mutate(
    p_chave = paste(dt_nascimento, mun_resid, sexo, sep = "_"),
    pseudo_id = digest(p_chave, algo = "md5")
  ) |>
  ungroup()

# identificando individuos que talvez estejam em meses diferentes

id_multt <- inss_ps |>
  group_by(pseudo_id) |>
  summarise(meses_p = n_distinct(competencia_concessao),
            regist = n()) |>
  filter(meses_p > 1)
# Visualizando a quantidade de individuos que estão em meses diferentes e a quantidade de vezes em que eles aparecem

p_dist_r <- id_multt %>%
  filter(meses_p > 1) %>%
  ggplot(aes(x = factor(meses_p))) +
  geom_bar(fill = "#2b5c8f", color = "black") +
  labs(title = "Distribuição de Indivíduos Recorrentes por Quantidade de Meses", x = "Quantidade de Meses Presente", y = "Número de Indivíduos (pseudo_id)") +
  theme_minimal()

recorrente <- inss_ps |>
  semi_join(id_multt |>
              filter(meses_p > 1), by = "pseudo_id")
# ----
p_sexo_r <- recorrente %>%
  count(sexo) %>%
  plot_ly(
    labels = ~ sexo,
    values = ~ n,
    type = "pie",
    textinfo = "label+percent",
    marker = list(colors = c("#2c7fb8", "#f03b20"))
  ) %>%
  layout(title = "Distribuição por sexo de indivíduos que receberam mais de um benefício no ano de 2025")
# p_sexo_r
top_mun_r <- recorrente %>%
  count(mun_resid, sort = TRUE) %>%
  slice_max(n, n = 15) %>%
  mutate(mun_resid = fct_reorder(mun_resid, n))

p_mun_r <- ggplot(top_mun_r, aes(x = n, y = mun_resid)) +
  geom_col(fill = "#756bb1") +
  labs(title = "Municípios com mais concessões de indivíduos mais de uma vez beneficiados no ano de 2025 (top 15)", x = "Nº de benefícios", y = NULL) +
  theme_minimal()
# p_mun_r

top_cid_r <- recorrente %>%
  count(cid_1, sort = TRUE) %>%
  slice_max(n, n = 15) %>%
  mutate(cid_1 = fct_reorder(cid_1, n))

p_cid_r <- ggplot(top_cid_r, aes(
  x = n,
  y = cid_1,
  text = paste0(cid_1, ": ", n)
)) +
  geom_col(fill = "#2c7fb8") +
  labs(title = "CIDs mais frequentes entre os beneficiados mais repetidos", x = "Nº de benefícios", y = NULL) +
  theme_minimal()
# p_cid_r

#----
# Verificando o banco apenas no ano de 2025, para beneficios iniciados no ano de 2025.


inss_ap_25 <- ap_inss_25 |>
  rowwise() |>
  mutate(
    p_chave = paste(dt_nascimento, mun_resid, sexo, sep = "_"),
    pseudo_id = digest(p_chave, algo = "md5")
  ) |>
  ungroup()

# identificando individuos que talvez estejam em meses diferentes

id_multt_ap <- inss_ap_25 |>
  group_by(pseudo_id) |>
  summarise(meses_p = n_distinct(competencia_concessao),
            regist = n()) |>
  filter(meses_p > 1)


p_dist_r1 <- id_multt_ap %>%
  filter(meses_p > 1) %>%
  ggplot(aes(x = factor(meses_p))) +
  geom_bar(fill = "#2b5c8f", color = "black") +
  labs(title = "Distribuição de Indivíduos Recorrentes por Quantidade de Meses", x = "Quantidade de Meses Presente", y = "Número de Indivíduos (pseudo_id)") +
  theme_minimal()

recorrente1 <- inss_ap_25 |>
  semi_join(id_multt_ap |>
              filter(meses_p > 1), by = "pseudo_id")

p_sexo_r1 <- recorrente1 %>%
  count(sexo) %>%
  plot_ly(
    labels = ~ sexo,
    values = ~ n,
    type = "pie",
    textinfo = "label+percent",
    marker = list(colors = c("#2c7fb8", "#f03b20"))
  ) %>%
  layout(title = "Distribuição por sexo de indivíduos que receberam mais de um benefício no ano de 2025")
# p_sexo_r
top_mun_r1 <- recorrente1 %>%
  count(mun_resid, sort = TRUE) %>%
  slice_max(n, n = 15) %>%
  mutate(mun_resid = fct_reorder(mun_resid, n))

p_mun_r1 <- ggplot(top_mun_r1, aes(x = n, y = mun_resid)) +
  geom_col(fill = "#756bb1") +
  labs(title = "Municípios com mais concessões de indivíduos mais de uma vez beneficiados no ano de 2025 (top 15)", x = "Nº de benefícios", y = NULL) +
  theme_minimal()
# p_mun_r

top_cid_r1 <- recorrente1 %>%
  count(cid_1, sort = TRUE) %>%
  slice_max(n, n = 15) %>%
  mutate(cid_1 = fct_reorder(cid_1, n))

p_cid_r1 <- ggplot(top_cid_r1, aes(
  x = n,
  y = cid_1,
  text = paste0(cid_1, ": ", n)
)) +
  geom_col(fill = "#2c7fb8") +
  labs(title = "CIDs mais frequentes entre os beneficiados mais repetidos", x = "Nº de benefícios", y = NULL) +
  theme_minimal()