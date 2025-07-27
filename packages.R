
required_packages <- c("readxl", "glmnet", "survival", "survminer", "gridExtra", "knitr")
install_if_missing <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
sapply(required_packages, install_if_missing)
