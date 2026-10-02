library(icenReg)
library(dplyr)
library(lubridate)
library(janitor)
library(arrow)
library(survival)
library(survminer)
# ----
cov_data_vec <- function(x) {
  # """
  # Converte vetor de strings em formato Date.
  # Aceita: datas em Excel, ISO (YYYY-MM-DD), europeu (DD/MM/YYYY), 
  # americano (MM/DD/YYYY), e com pontos/traços.
  # 
  # Totalmente vetorizada - compatível com dplyr::across()
  # Sem avisos de coerção.
  # """
  
  result <- rep(NA_Date_, length(x))
  valid_idx <- !is.na(x) & x != ""
  
  if (!any(valid_idx)) return(result)
  
  x_trimmed <- trimws(x[valid_idx])
  
  # ---- PASSO 1: Identifica datas numéricas (Excel) ----
  is_numeric <- grepl("^[0-9]{4,6}$", x_trimmed)
  numeric_idx <- which(valid_idx)[is_numeric]
  
  if (any(is_numeric)) {
    num_vals <- suppressWarnings(as.numeric(x_trimmed[is_numeric]))
    valid_excel <- !is.na(num_vals) & num_vals < 60000
    
    result[numeric_idx[valid_excel]] <- 
      janitor::excel_numeric_to_date(num_vals[valid_excel])
  }
  
  # ---- PASSO 2: Processa datas textuais ----
  text_idx <- which(valid_idx)[!is_numeric]
  if (length(text_idx) > 0) {
    text_dates <- suppressWarnings(
      as.Date(
        lubridate::parse_date_time(
          x_trimmed[!is_numeric],
          orders = c("ymd", "dmy", "mdy", "Ymd"),
          quiet = TRUE
        )
      )
    )
    result[text_idx] <- text_dates
  }
  
  return(result)
}

inss_pa <- read_parquet("inss_pa.parquet") |>
  clean_names() |>
  mutate(
    across(starts_with("dt_"), ~ cov_data_vec(.x))
  ) 

datas <- inss_pa |>
  filter(dt_dib >= "2025-01-01" & dt_dib <= "2025-12-31") |>
  mutate( 
    espera = round(time_length(interval(dt_dib, dt_ddb), "month")),
    status = if_else(is.na(dt_ddb), 0, 1),
    duracao = round(time_length(interval(dt_dib, dt_dcb),"month"))
  ) |>
  filter(
    espera >= 0, duracao >= 0
  )
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
# ----

# fit_0 <- ic_np(cbind(left,right) ~ 1, data = datas)
# summary(fit_0)

# evento: tempo até a concessão do benefício
fit_surv <- survfit(Surv(espera, status) ~ sexo, data = datas)
summary(fit_surv)

# evento: tempo de duração do período de afastamento 
# sexo
# entre forma de filiação (desempregado, empregado,doméstico, autonomo, facultativo(contribui, desempregado), segurado especial(trab. rural e afins) e trab. avulso (sem vinculo, pretador de serviço com sindicato))
# ramo de atividade
fit_surv_s <- survfit(Surv(duracao) ~ ramo_atividade, data = datas)
summary(fit_surv_s)


plot(fit_surv, xlab = "Tempo (dias)", ylab = "Sobrevivência")

ggsurvplot(fit_surv,
            datas,
            fun = "pct", #Exibe a porcentagem
            palette = c("darkblue","darkred"),
            conf.int = T,
            cumcensor = T,
            tables.height = 0.25,
            xlab = "Tempo (em meses)",
            ylab = "Sobrevivência Livre do Evento",
            title = "Tempo até a concessão do benefício"
)

ggsurvplot(fit_surv_s$strata[-2],
           datas,
           fun = "pct", #Exibe a porcentagem
           conf.int = T,
           cumcensor = T,
           tables.height = 0.25,
           xlab = "Tempo (em meses)",
           ylab = "Sobrevivência Livre do Evento",
           title = "Tempo de duração de afastamento"
)
