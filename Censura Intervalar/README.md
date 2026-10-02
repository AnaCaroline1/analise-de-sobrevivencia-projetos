# 📚 Estimador de Turnbull com ggplot - Guia Completo

## 📖 O Que Você Recebeu

Este pacote contém **tudo** o que você precisa para plotar o Estimador de Turnbull com ggplot, comparar grupos e adicionar intervalos de confiança.

---

## 📁 Arquivos Inclusos

### 1️⃣ **COMECE AQUI** 
- **`ANTES_DEPOIS_VISUALIZACAO.md`** ← **LEIA PRIMEIRO**
  - Explica o erro comum (gráfico em 100%)
  - Mostra a diferença entre código errado e correto
  - Tem diagramas visuais do problema e solução
  - ⏱️ Leitura: 5-10 minutos

### 2️⃣ **CÓDIGOS PRONTOS PARA USAR**

- **`turnbull_plot_ggplot.R`** - Exemplo completo e genérico
  - Cria dados fictícios de exemplo
  - Aplica Turnbull em dois grupos
  - Gera 3 gráficos diferentes com ggplot
  - Mostra tabelas de resultados
  - ✅ **Copie e adapte para seus dados**

- **`turnbull_exemplo_cancer_mama.R`** - Exemplo realista (seu caso)
  - Simula dados de tempo até recorrência de câncer de mama
  - Dois grupos de tratamento (Radioterapia vs Quimio+Rad)
  - Gera gráficos publicáveis
  - Inclui análise estatística descritiva
  - ✅ **Próximo à estrutura do seu `c_mama`**

- **`CHEAT_SHEET_RAPIDO.R`** - Códigos para copiar e colar
  - 4 opções de gráficos (da mais rápida à mais bonita)
  - Cores recomendadas
  - Customizações comuns
  - ✅ **Para quando você quer resultados rápidos**

### 3️⃣ **DOCUMENTAÇÃO**

- **`GUIA_TURNBULL_GGPLOT.md`** - Guia completo e técnico
  - Estrutura de dados necessária
  - Passo a passo da plotagem
  - 3 opções de visualização
  - Explicação de cada elemento do gráfico
  - Quando usar cada tipo de gráfico

- **`TROUBLESHOOTING_RAPIDO.md`** - Resolução de problemas
  - 10+ problemas comuns e soluções
  - Causas raiz de cada erro
  - Código correto vs errado
  - Checklist final
  - ✅ **Consulte quando algo der errado**

### 4️⃣ **ARQUIVO ESSENCIAL**

- **`turnbull_estimator.R`** - Seu estimador de Turnbull
  - Funções: `turnbull_em()`, `turnbull_table()`, `turnbull_boot_ci()`
  - Use com: `source("turnbull_estimator.R")`

---

## 🚀 INÍCIO RÁPIDO (3 passos)

### Passo 1: Entender o Problema (5 min)
```bash
Leia: ANTES_DEPOIS_VISUALIZACAO.md
```

### Passo 2: Ver um Exemplo (10 min)
```r
# No RStudio, abra e rode:
source("turnbull_plot_ggplot.R")
# Vai gerar 3 gráficos automaticamente
```

### Passo 3: Adaptar para Seus Dados (30 min)
```r
# Copie e adapte turnbull_exemplo_cancer_mama.R
# Substitua seus dados
# Rode!
```

---

## 📊 Qual Gráfico Usar?

| Situação | Use | Arquivo |
|----------|-----|---------|
| Quero entender o problema | `ANTES_DEPOIS_VISUALIZACAO.md` | - |
| Quero código rápido | `CHEAT_SHEET_RAPIDO.R` | Opção 1 |
| Quero gráfico bonito | `CHEAT_SHEET_RAPIDO.R` | Opção 2 |
| Quero gráficos separados | `CHEAT_SHEET_RAPIDO.R` | Opção 3 |
| Quero exemplo completo | `turnbull_plot_ggplot.R` | - |
| Quero case de câncer mama | `turnbull_exemplo_cancer_mama.R` | - |
| Tenho um problema | `TROUBLESHOOTING_RAPIDO.md` | - |
| Quero aprender detalhes | `GUIA_TURNBULL_GGPLOT.md` | - |

---

## ✅ Estrutura Correta (O Que Você Estava Perdendo)

### ❌ Errado (seu problema original)
```r
fit_1 <- turnbull_em(L, U)
# ❌ Tentava plotar fit_1 direto
# ❌ Resultado: gráfico em 100%
```

### ✅ Correto (agora você sabe!)
```r
fit_1 <- turnbull_em(L, U)
tab_1 <- turnbull_table(fit_1)  # ← CRUCIAL!
# ✅ Plota tab_1 com ggplot
# ✅ Resultado: curva de sobrevida com IC
```

---

## 📋 Checklist: Você tem tudo?

- ✅ Arquivo `turnbull_estimator.R` (seu código original)
- ✅ Script `turnbull_plot_ggplot.R` (exemplo genérico)
- ✅ Script `turnbull_exemplo_cancer_mama.R` (exemplo realista)
- ✅ Script `CHEAT_SHEET_RAPIDO.R` (códigos prontos)
- ✅ Documento `ANTES_DEPOIS_VISUALIZACAO.md` (entender o erro)
- ✅ Documento `GUIA_TURNBULL_GGPLOT.md` (aprender)
- ✅ Documento `TROUBLESHOOTING_RAPIDO.md` (resolver problemas)

---

## 🎯 Roadmap: Por Onde Começo?

### Se você quer **entender o problema agora**:
1. Leia `ANTES_DEPOIS_VISUALIZACAO.md` (10 min)
2. Olhe o diagrama de fluxo correto
3. Veja a comparação "❌ vs ✅"

### Se você quer **código funcionando agora**:
1. Abra `CHEAT_SHEET_RAPIDO.R`
2. Copie a "OPÇÃO 2: GRÁFICO BONITO"
3. Adapte para seus dados (5 min)
4. Rode!

### Se você quer **aprender direito**:
1. Leia `ANTES_DEPOIS_VISUALIZACAO.md` (10 min)
2. Rode `turnbull_exemplo_cancer_mama.R` (5 min)
3. Leia `GUIA_TURNBULL_GGPLOT.md` (20 min)
4. Adapte para seus dados

### Se deu erro:
1. Leia `TROUBLESHOOTING_RAPIDO.md`
2. Procure seu erro específico
3. Copie a solução fornecida

---

## 🔑 Conceitos-Chave

### A Função Crítica
```r
turnbull_table(fit, conf.level = 0.95)
```
Retorna:
- `survival` → A curva de sobrevida que você quer plotar
- `lower95`, `upper95` → Intervalo de confiança 95%
- `n.risk` → Número em risco
- `n.event` → Número de eventos esperados

### O Gráfico Correto
```r
ggplot(tab, aes(x = p, y = survival)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.2) +
  geom_step(direction = "hv") +  # ← Chave: geom_step, não geom_line!
  scale_y_continuous(labels = scales::percent)
```

### Passos Sempre Nessa Ordem
1. `L`, `U` (dados brutos)
2. `turnbull_em(L, U)` (ajuste do modelo)
3. `turnbull_table(fit)` (extrair sobrevida + IC) ← **CRUCIAL**
4. `ggplot()` (plotar)

---

## 📞 Dúvidas Frequentes

### "Por que meu gráfico fica em 100%?"
→ Você não está usando `turnbull_table()`. Leia `ANTES_DEPOIS_VISUALIZACAO.md`.

### "Como comparar dois grupos?"
→ Combine as duas tabelas de `turnbull_table()` em um só `data.frame`. Veja `turnbull_exemplo_cancer_mama.R`.

### "Posso usar bootstrap para IC?"
→ Sim! Use `turnbull_boot_ci()`. Veja `CHEAT_SHEET_RAPIDO.R` - OPÇÃO: USO IC POR BOOTSTRAP.

### "Qual é a diferença entre `turnbull_em()` e `turnbull_table()`?"
→ Leia a tabela em `ANTES_DEPOIS_VISUALIZACAO.md` - seção "A Diferença Chave".

### "Como salvar o gráfico?"
→ Use `ggsave()`. Exemplos em `CHEAT_SHEET_RAPIDO.R`.

### "Como personalizar cores?"
→ Veja `CHEAT_SHEET_RAPIDO.R` - seção "CORES BONITAS".

---

## 🎓 Leitura Recomendada (por ordem)

1. **5 min** → `ANTES_DEPOIS_VISUALIZACAO.md` (entender problema)
2. **10 min** → Rode `turnbull_exemplo_cancer_mama.R` (ver funcionando)
3. **20 min** → `GUIA_TURNBULL_GGPLOT.md` (aprender detalhes)
4. **Conforme necessário** → `TROUBLESHOOTING_RAPIDO.md` (resolver erros)

---

## 💡 Dica de Ouro

**Comece com o código mais simples** (`CHEAT_SHEET_RAPIDO.R` - OPÇÃO 1), entenda como funciona, depois adicione complexity (IC, cores, customizações).

```r
# Simples
ggplot(tab, aes(x = p, y = survival)) +
  geom_step(direction = "hv") +
  scale_y_continuous(labels = scales::percent)

# Meio
# ... + geom_ribbon() para IC
# ... + geom_point() para pontos

# Complexo
# ... + colors, facets, temas, legenda posicionada, etc
```

---

## ✨ Próximos Passos

Após usar este material:

1. ✅ Adapte para seus dados reais (`c_mama`)
2. ✅ Gere gráficos para publicação
3. ✅ Salve em PNG (300 dpi) e PDF
4. ✅ Crie tabela de resultados em CSV
5. ✅ Escreva a seção de Resultados do seu artigo/tese

---

## 📝 Resumo Rápido

| O que fazer | Onde está |
|------------|----------|
| Entender o erro | `ANTES_DEPOIS_VISUALIZACAO.md` |
| Copiar código | `CHEAT_SHEET_RAPIDO.R` |
| Ver exemplo | `turnbull_exemplo_cancer_mama.R` |
| Aprender detalhes | `GUIA_TURNBULL_GGPLOT.md` |
| Resolver problema | `TROUBLESHOOTING_RAPIDO.md` |

---

## 🎉 Você está pronto!

Abra `ANTES_DEPOIS_VISUALIZACAO.md` e comece agora! 📊✨

---

**Versão:** 1.0  
**Data:** Setembro 2026  
**Baseado em:** Estimador de Turnbull (NPMLE) para censura intervalar + ggplot2
