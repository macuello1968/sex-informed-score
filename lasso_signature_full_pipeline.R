# ============================================================================
# Script: lasso_signature_full_pipeline.R
# Title: Full Pipeline for LASSO-Based Prognostic Score Construction and Validation
# Author: Mauricio Cuello Fredes
# Description:
#   Complete pipeline to derive and validate a prognostic gene expression
#   signature using LASSO Cox regression on Overall Survival (OS) and
#   Progression-Free Interval (PFI).
#
#   Steps:
#     1) Load z-scored gene expression and survival data
#     2) Run LASSO Cox for OS and PFI
#     3) Keep genes selected in BOTH models
#     4) Build a 10-gene score using the mean coefficient (OS/PFI)
#     5) Dichotomize patients (High/Low) by median
#     6) Produce Kaplan–Meier curves for OS & PFI
#     7) Compare 10-gene vs 9-gene score (excluding CGB8)
#     8) Print the explicit score formula with coefficients
#     9) Save plots to /outputs
#
# Input:
#   - Excel file with:
#       * 30 candidate genes (Z-scores)
#       * OS, OS_months, PFI, PFI_months
#   File name expected: LASSO_definitive_dataset.xlsx
#
# Where to put the dataset (GitHub-friendly):
#   <repo_root>/data/LASSO_definitive_dataset.xlsx
#
# Outputs:
#   - Printed list of selected genes
#   - Printed score formula with gene-specific coefficients
#   - KM plots (OS/PFI) for 10-gene and 9-gene scores
#   - PDFs saved under ./outputs/
#
# Reproducibility:
#   R >= 4.1.0
#   Packages: readxl, glmnet, survival, survminer, gridExtra
#
# Usage:
#   - Place the dataset under ./data/LASSO_definitive_dataset.xlsx
#   - Run: source("lasso_signature_full_pipeline.R")
# ============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(glmnet)
  library(survival)
  library(survminer)
  library(gridExtra)
})
set.seed(12345)

# ==== 1) Path Config ====
if (Sys.info()[["user"]] == "mauriciocuello") {
  data_path <- "/Users/mauriciocuello/Desktop/LASSO_definitive_dataset.xlsx"
} else {
  data_path <- file.path("data", "LASSO_definitive_dataset.xlsx")
}

out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

if (!file.exists(data_path)) {
  stop(paste0("Input file not found: ", data_path))
}

# ==== 2) Load Data ====
df <- read_excel(data_path)

# ==== 3) Define Genes ====
genes_z <- c("ZABCA13", "ZADAM12", "ZADAM15", "ZADAMTS5", "ZANO4", "ZAP1S1", "ZAPBB3", "ZAPC2", "ZARID1A",
             "ZB4GALNT1", "ZCD276", "ZCGB8", "ZCOL11A1", "ZCOL14A1", "ZCOL25A1", "ZEIF3J", "ZFAT1",
             "ZIGF2BP2", "ZIKBKE", "ZIRX3", "ZPHOX2A", "ZPIK3CA", "ZPPM1L", "ZPPP1R13B", "ZPYGM",
             "ZRASIP1", "ZRNF43", "ZSAMD9", "ZTEK", "ZTERT")

# ==== 4) LASSO Cox Models ====
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

# ==== 5) Final Genes & Score ====
genes_final <- intersect(genes_os, genes_pfi)
coefs_os_vec <- as.vector(coefs_os[genes_final, ])
coefs_pfi_vec <- as.vector(coefs_pfi[genes_final, ])
coef_mean <- (coefs_os_vec + coefs_pfi_vec) / 2
names(coef_mean) <- genes_final

expr_all <- as.matrix(df[, genes_final])
df$score_10gene <- as.numeric(expr_all %*% coef_mean)
cutoff <- median(df$score_10gene, na.rm = TRUE)
df$score_group <- ifelse(df$score_10gene >= cutoff, "High", "Low")

# ==== 6) Print Genes and Formula ====
cat("✅ Genes retained:
")
print(genes_final)
cat("
🧮 score_10gene[i] = sum(beta_j * Z_ij)

")
cat("score_10gene = ", paste(sprintf("(%.4f × %s)", coef_mean, names(coef_mean)), collapse = " + "), "

")

# ==== 7) KM Curves (10-gene) ====
df_km_os <- df[complete.cases(df$OS_months, df$OS, df$score_group) & df$OS_months > 0, ]
fit_os <- survfit(Surv(OS_months, OS) ~ score_group, data = df_km_os)
plot_os <- ggsurvplot(fit_os, data = df_km_os, pval = TRUE, risk.table = TRUE)
print(plot_os)
ggsave(file.path(out_dir, "KM_OS_10gene.pdf"), plot = plot_os$plot)

df_km_pfi <- df[complete.cases(df$PFI_months, df$PFI, df$score_group) & df$PFI_months > 0, ]
fit_pfi <- survfit(Surv(PFI_months, PFI) ~ score_group, data = df_km_pfi)
plot_pfi <- ggsurvplot(fit_pfi, data = df_km_pfi, pval = TRUE, risk.table = TRUE)
print(plot_pfi)
ggsave(file.path(out_dir, "KM_PFI_10gene.pdf"), plot = plot_pfi$plot)

# ==== 8) KM Comparison: 10-gene vs 9-gene (excl. CGB8) ====
genes_9 <- setdiff(names(coef_mean), "ZCGB8")
df$score_9gene <- rowMeans(df[, genes_9], na.rm = TRUE)
df$group_10 <- ifelse(df$score_10gene >= median(df$score_10gene, na.rm = TRUE), "High", "Low")
df$group_9 <- ifelse(df$score_9gene >= median(df$score_9gene, na.rm = TRUE), "High", "Low")

df_km <- df[complete.cases(df$OS_months, df$OS, df$PFI_months, df$PFI, df$group_10, df$group_9) &
              df$OS_months > 0 & df$PFI_months > 0, ]
fit_os_10 <- survfit(Surv(OS_months, OS) ~ group_10, data = df_km)
fit_os_9 <- survfit(Surv(OS_months, OS) ~ group_9, data = df_km)
fit_pfi_10 <- survfit(Surv(PFI_months, PFI) ~ group_10, data = df_km)
fit_pfi_9 <- survfit(Surv(PFI_months, PFI) ~ group_9, data = df_km)

p1 <- ggsurvplot(fit_os_10, data = df_km, pval = TRUE)
p2 <- ggsurvplot(fit_os_9, data = df_km, pval = TRUE)
pdf(file.path(out_dir, "KM_OS_comparison.pdf"))
grid.arrange(p1$plot, p2$plot, ncol = 2)
dev.off()

p3 <- ggsurvplot(fit_pfi_10, data = df_km, pval = TRUE)
p4 <- ggsurvplot(fit_pfi_9, data = df_km, pval = TRUE)
pdf(file.path(out_dir, "KM_PFI_comparison.pdf"))
grid.arrange(p3$plot, p4$plot, ncol = 2)
dev.off()

# ==== 9) Session Info ====
cat("

--- Session info ---
")
print(sessionInfo())
