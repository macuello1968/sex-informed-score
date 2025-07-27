
# === 0. LIBRERÍAS NECESARIAS ===
library(readxl)
library(glmnet)
library(survival)
library(survminer)
library(gridExtra)

# === 1. CARGAR BASE ===
ruta <- "/Users/mauriciocuello/Desktop/CCR_score/LASSO_definitive_ethnicity_simplified.xlsx"
df <- read_excel(ruta)

# === 2. DEFINIR LOS 30 GENES Z-SCORE ===
genes_z <- c("ZABCA13", "ZADAM12", "ZADAM15", "ZADAMTS5", "ZANO4", "ZAP1S1", "ZAPBB3", "ZAPC2", "ZARID1A",
             "ZB4GALNT1", "ZCD276", "ZCGB8", "ZCOL11A1", "ZCOL14A1", "ZCOL25A1", "ZEIF3J", "ZFAT1",
             "ZIGF2BP2", "ZIKBKE", "ZIRX3", "ZPHOX2A", "ZPIK3CA", "ZPPM1L", "ZPPP1R13B", "ZPYGM",
             "ZRASIP1", "ZRNF43", "ZSAMD9", "ZTEK", "ZTERT")

# === 3-10: Filtrar, ajustar modelos LASSO, calcular scores ===
df_os <- df[complete.cases(df$OS_months, df$OS) & df$OS_months > 0, ]
expr_os <- as.matrix(df_os[, genes_z])
cox_os <- cv.glmnet(expr_os, Surv(df_os$OS_months, df_os$OS), family = "cox", alpha = 1)
coefs_os <- coef(cox_os, s = "lambda.min")
genes_os <- rownames(coefs_os)[which(coefs_os != 0)][-1]

df_pfi <- df[complete.cases(df$PFI_months, df$PFI) & df$PFI_months > 0, ]
expr_pfi <- as.matrix(df_pfi[, genes_z])
cox_pfi <- cv.glmnet(expr_pfi, Surv(df_pfi$PFI_months, df_pfi$PFI), family = "cox", alpha = 1)
coefs_pfi <- coef(cox_pfi, s = "lambda.min")
genes_pfi <- rownames(coefs_pfi)[which(coefs_pfi != 0)][-1]

genes_final <- intersect(genes_os, genes_pfi)
coefs_os_vec <- as.vector(coefs_os[genes_final, ])
coefs_pfi_vec <- as.vector(coefs_pfi[genes_final, ])
coef_mean <- (coefs_os_vec + coefs_pfi_vec) / 2
names(coef_mean) <- genes_final
expr_all <- as.matrix(df[, genes_final])
df$score_10gene <- as.numeric(expr_all %*% coef_mean)
cutoff <- median(df$score_10gene, na.rm = TRUE)
df$score_group <- ifelse(df$score_10gene >= cutoff, "High", "Low")

# === 11–12: Kaplan–Meier ===
df_km_os <- df[complete.cases(df$OS_months, df$OS, df$score_group) & df$OS_months > 0, ]
fit_os <- survfit(Surv(OS_months, OS) ~ score_group, data = df_km_os)
plot_os <- ggsurvplot(fit_os, data = df_km_os, pval = TRUE, risk.table = TRUE,
                      title = "Overall Survival según score LASSO (30 → 10 genes)",
                      xlab = "Meses", ylab = "Probabilidad de sobrevida")
ggsave("~/Desktop/OS_10gene.pdf", plot = plot_os$plot, width = 6, height = 5)

df_km_pfi <- df[complete.cases(df$PFI_months, df$PFI, df$score_group) & df$PFI_months > 0, ]
fit_pfi <- survfit(Surv(PFI_months, PFI) ~ score_group, data = df_km_pfi)
plot_pfi <- ggsurvplot(fit_pfi, data = df_km_pfi, pval = TRUE, risk.table = TRUE,
                       title = "Progression-Free Interval según score LASSO (30 → 10 genes)",
                       xlab = "Meses", ylab = "Probabilidad libre de progresión")
ggsave("~/Desktop/PFI_10gene.pdf", plot = plot_pfi$plot, width = 6, height = 5)
