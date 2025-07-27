# Sex-Informed Prognostic Score Based on LASSO Cox Models

This repository contains the full R pipeline and dataset used to generate and validate a sex-informed 10-gene prognostic signature using LASSO Cox regression. The model was developed using z-score normalized gene expression data and survival outcomes (Overall Survival and Progression-Free Interval), and is applicable across multiple cancer types.

## 🔍 Summary

- Selects prognostic genes using LASSO Cox models for OS and PFI
- Identifies genes shared by both models
- Builds a composite prognostic score (10-gene and 9-gene versions)
- Stratifies patients into High vs. Low risk groups
- Generates Kaplan–Meier survival curves
- Outputs the final score formulas with gene-specific coefficients

## 📁 Folder Structure

```
sex-informed-score/
├── data/
│   └── LASSO_definitive_dataset.xlsx       # Input dataset (not public)
├── lasso_signature_full_pipeline.R         # Main R script
├── outputs/                                # Generated plots
├── README.md                               # Project documentation
├── LICENSE                                 # MIT license
├── LASSO_Readme_Documentation.docx         # Supplementary documentation
```

## 💻 Requirements

- R version ≥ 4.1.0
- Install required R packages:

```r
install.packages(c("readxl", "glmnet", "survival", "survminer", "gridExtra"))
```

## 🚀 How to Run

1. Clone the repository or download the files manually.
2. Place the dataset in the `data/` folder as:
   ```
   data/LASSO_definitive_dataset.xlsx
   ```
   Alternatively, edit the `data_path` variable in the script to reflect your absolute path.

3. Open R or RStudio and run:
   ```r
   source("lasso_signature_full_pipeline.R")
   ```

4. The script will print the selected genes, the coefficients, and save Kaplan–Meier plots in the `outputs/` folder.

## 📊 Output

- Printed formulas for:
  - 10-gene score
  - 9-gene score (excluding CGB8)
- Kaplan–Meier plots:
  - OS_10gene.pdf
  - PFI_10gene.pdf
  - OS_9gene.pdf
  - PFI_9gene.pdf
  - KM_OS_comparison.pdf
  - KM_PFI_comparison.pdf

## 📌 Reproducibility

- The script uses `set.seed(12345)` for reproducible LASSO selection
- Final session details printed using `sessionInfo()`

## 📄 License

This repository is released under the MIT License. It is intended for academic and non-commercial use. Please cite the corresponding manuscript when using this code, dataset, or methodology.

## 📬 Contact

**Mauricio Cuello Fredes**  
Pontificia Universidad Católica de Chile  
✉️ mcuello@uc.cl
