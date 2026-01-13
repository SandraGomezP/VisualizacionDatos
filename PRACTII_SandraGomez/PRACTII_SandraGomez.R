# ==============================================================================
# 1. CONFIGURACIÓN Y LIBRERÍAS
# ==============================================================================
library(tidyverse)
library(janitor)
library(survival)
library(survminer)
library(plotly)
library(data.table)

# Configurar directorio (ajusta según tu PC)
setwd("D:/SANDRA/UOC/Visualización/PRACT")

# ==============================================================================
# 2. CARGA DE DATOS (CLÍNICOS, RPPA Y MUTACIONES)
# ==============================================================================
# Datos Clínicos
df_raw <- read_tsv("brca_tcga_pan_can_atlas_2018_clinical_data.tsv", comment = "#")

# Datos de Proteínas (RPPA)
rppa_raw <- fread("data_rppa.txt", data.table = FALSE)

# Datos de Mutaciones
mutations <- fread("data_mutations.txt", data.table = FALSE)

# ==============================================================================
# 3. PROCESAMIENTO Y LIMPIEZA INTEGRADA
# ==============================================================================

# A. Limpieza de Clínicos
df_clean <- df_raw %>%
  clean_names() %>%
  mutate(
    status_binary = ifelse(overall_survival_status == "1:DECEASED", 1, 0),
    aneuploidy_level = ifelse(aneuploidy_score > median(aneuploidy_score, na.rm = TRUE), 
                              "High", "Low")
  ) %>%
  filter(subtype %in% c("BRCA_LumA", "BRCA_LumB", "BRCA_Basal", "BRCA_Her2"))

# B. Procesamiento de Proteínas (Transposición y Nombres Únicos)
proteinas <- rppa_raw[, 1]
proteinas[is.na(proteinas) | proteinas == ""] <- "Unknown_Protein"
proteinas_unicas <- make.unique(proteinas)

rppa_t <- as.data.frame(t(rppa_raw[, -1]))
colnames(rppa_t) <- proteinas_unicas
rppa_final <- rppa_t %>% as_tibble(rownames = "sample_id")

# C. Identificación de Mutaciones en PIK3CA
pik3ca_mutated_samples <- mutations %>%
  filter(Hugo_Symbol == "PIK3CA") %>%
  mutate(sample_id_short = substr(Tumor_Sample_Barcode, 1, 15)) %>%
  pull(sample_id_short) %>%
  unique()

# ==============================================================================
# 4. UNIÓN FINAL DE DATASETS (MASTER TABLE)
# ==============================================================================
df_completo <- df_clean %>%
  inner_join(rppa_final, by = "sample_id") %>%
  mutate(
    pik3ca_status = ifelse(sample_id %in% pik3ca_mutated_samples, "Mutado", "Normal"),
    trastuzumab_group = ifelse(grepl("Trastuzumab", neoadjuvant_therapy_type_administered_prior_to_resection_text, ignore.case = TRUE), 
                               "Trastuzumab", "Otros/No"),
    dfs_event = ifelse(disease_free_status == "1:Recurred/Progressed", 1, 0)
  )

# ==============================================================================
# 5. GENERACIÓN DE RESULTADOS (Q1, Q2, Q3)
# ==============================================================================

### [Q1] Supervivencia Global por Aneuploidía (Subtipo Luminal B)
df_luminalB <- df_completo %>% filter(subtype == "BRCA_LumB")
fit_q1 <- survfit(Surv(overall_survival_months, status_binary) ~ aneuploidy_level, data = df_luminalB)

p1 <- ggsurvplot(fit_q1, data = df_luminalB, risk.table = TRUE, pval = TRUE,
                 conf.int = TRUE, palette = c("#E7B800", "#2E9FDF"),
                 title = "Q1: Supervivencia Global vs Aneuploidía (Luminal B)",
                 xlab = "Meses")
p1

### [Q2] Proteína HER2 vs Metástasis (Subtipo Basal) - INTERACTIVO
p2_static <- df_completo %>%
  filter(subtype == "BRCA_Basal") %>%
  ggplot(aes(x = american_joint_committee_on_cancer_metastasis_stage_code, 
             y = `ERBB2|HER2_pY1248`, 
             fill = american_joint_committee_on_cancer_metastasis_stage_code,
             text = paste("Paciente:", patient_id))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  theme_minimal() +
  labs(title = "Q2: Niveles de HER2 Fosforilada en Subtipo Basal",
       x = "Estadio de Metástasis", y = "Z-score Proteína", fill = "Estadio")

ggplotly(p2_static, tooltip = "text")

### [Q3] Respuesta a Terapia (DFS): Trastuzumab e Interacción Genética
df_q3 <- df_completo %>% 
  filter(!is.na(disease_free_months)) %>%
  mutate(grupo_combinado = paste(trastuzumab_group, "-", pik3ca_status))

fit_q3 <- survfit(Surv(disease_free_months, dfs_event) ~ grupo_combinado, data = df_q3)

p3 <- ggsurvplot(fit_q3, data = df_q3, pval = TRUE, risk.table = TRUE, palette = "Set1",
                 title = "Q3: Supervivencia Libre de Enfermedad (DFS)",
                 subtitle = "Interacción Trastuzumab - Mutación PIK3CA",
                 xlab = "Meses", legend.title = "Grupo")
p3

# ==============================================================================
# 6. ANÁLISIS EXTRA: RELACIÓN EDAD VS ANEUPLOIDÍA
# ==============================================================================
ggplot(df_completo, aes(x = diagnosis_age, y = aneuploidy_score, color = subtype)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~subtype) +
  theme_bw() +
  labs(title = "Relación Edad vs Inestabilidad Genómica",
       x = "Edad al Diagnóstico", y = "Aneuploidy Score")
