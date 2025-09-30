# Polygenic Score and Disease Prevalence Analysis

## Abstract
This project investigates the correlation between polygenic scores (PGS) and
the prevalence of multiple complex diseases such as type 2 diabetes mellitus,
chronic kidney disease, breast cancer, and others, across diverse human populations.  
Please refer to the file **report_country_pgs.pdf** for more details.

---

## Project Structure

### 📂 results
- `logistic_model_summary.csv` – Summary results of logistic regression models (for disease).   
- `fig_by_cause_logit/` – Figures generated from logistic regression, grouped by disease.  
- `fig_quant_traits/` – Figures for quantitative trait analyses.

### 📂 scores
- `PGSxxxxx_scores.txt` – Polygenic score files for different traits/diseases. Each file contains the computed scores with pgs_cal.

### 📂 codes
- `pgs_disease_logistic.R` – R script for logistic regression modeling on PGS vs. disease.  
- `lin_qt.R` – R script for linear models on PGS vs. quantitative traits.  
- `multiscore.sh` – Shell script for running multiple PGS analyses in batch.  
- `country_pgs.Rproj` – R project configuration file for convenient project management.

### 📂 data
- `cause_pgs_map.csv` – Mapping between diseases/traits and their associated PGS Catalog ID.
- `points_pgs_vs_bmi_female.csv` – Metadata for female BMI.
- `points_pgs_vs_bmi_male.csv` – Metadata for male BMI.
- `points_pgs_vs_height_female.csv` – Metadata for female Height
- `points_pgs_vs_height_male.csv` – Metadata for male Height.

the above 4 added mean PGS for populations as new coloum.

- `population_bmi.csv` – Country-level BMI data.  
- `population_height.csv` – Country-level height data.  
- `pop_super_metadata.csv` – Extracted metadata for samples in 1000 Genomes Project
- `IHME-GBD_2021_DATA-cff1c199-1.csv` – Global Burden of Disease (GBD) 2021 dataset used as reference.  
- `integrated_call_samples_v3.20130502.ALL.panel` – 1000 Genomes Project integrated sample panel information - from the website.  
- `20131219.populations.tsv` – Population and Country labels of the 1000 Genomes Project - from the website.

### 📄 report_country_pgs.pdf
Main project report with methodology, results, and interpretation.

---
