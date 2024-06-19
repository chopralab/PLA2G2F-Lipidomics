# README

# Overview
## This github repo contains code from the Comprehesive Lipidomic Analysis Workflow (CLAW), developed in the Chopra Laboratory at Purdue University. This workflow was used to investage the role of Pla2g2f in dentate gyrus homeostasis, neuronal health, and cognitive resilience. Pla2g2f was virally ablated in dentate granule cells (DGC), and the entire region was dissected and processed for lipidomics analysis using UHPLC-MSMS. All identified lipid species were exported into the excel file (.xls) in the data folder and used for downstream analysis, including differential expression and pathway enrichment. 

# Data
## All identified lipid species were exported into an Excel file (.xls) located in the data folder. This data serves as the basis for downstream analyses, including differential expression and pathway enrichment studies.

# Analysis Scripts
## The repository includes R scripts for lipidomic analysis, facilitating the processing and interpretation of the data. The three steps in analysis include data processing, EdgeR analysis, and data visualization.


# 1 Data Processing

# This code accomplishes 3 things
## 1. filter the raw data for grade A & B lipids. 
## 2. sum ion itensities for lipids with identical formulas (misidentified as two peaks in UHPLC)
## 3. create tables of the number of lipid species detected and total ion intensity detected in each lipid class for overview of data

# Lipid grading system used:
## A: All parent and acyl chains detected
## B: Parent and at least one acyl chain detected
## C: no parent but both acyl chains detected
## D: only single acyl chain detected
## Grade C & D were poorly matched and not advised to include in the analysis.

# This code outputs 6 files:
## 1. sum_identical_lipids.csv contains all 1493 unique lipids (grade A-D)
## 2. sum_identical_lipids_AB_only.csv contains 736 unique lipids (grade A&B only)
## 3. all_lipids_detected.csv contains the total number of lipid species detected in each lipid class (grade A-D)
## 4. grade_AB_lipids_detected.csv contains the total number of lipid species detected in each lipid class (grade A&B only)
## 5. sum_ion_intensity_per_class_all.csv contains the total signal (ion intensity) detected in each lipid class (grade A-D)
## 6. sum_ion_intensity_per_class_grade_AB.csv contains the total signal (ion intensity) detected in each lipid class (grade A&B only)

# The second output file (sum_identical_lipids_AB_only.csv) was used for downstream differential expression analysis. 

# 2 EdgeR Analysis
# This code performs EdgeR analysis for pre-processed data.
## EdgeR is a form of differential expression analysis which uses a generalized linear model (GLM). 
## The GLM uses a negative binomial distribution to model lipid ion count data (or mRNA count data) accounting for the technical and biological variability. 

## Pre-processing is required to ensure there is no duplicate lipid entries in the input data. 

# This code generates one output file: Pla2g2fko_vs_WT_gradeAB.csv
## Columns in this output file include:
### lipid: lipid name in this general format: Headgroup(sn1/sn2/sn3)
### type: lipid type (abbreviations can be found in the supplement materials: ________), 
### mean1 & mean2: mean ion intensity for the two groups (Group1: WT, Group2: Pla2g2f cKO)
### logFC: Log2 transformed fold change of Group2 (Pla2g2f cKO) over Group1 (WT)
### logCPM: Log2 transformed average expression of the lipid species across all samples
### LR: likelihood ratio test statistic
### PValue: nominal p-value derived from qualsi-likelihood ratio test without multiple testing correction
### FDR: false discovery rate calculated using the Benjamini–Hochberg method
### Z1_sum, Z2_sum...K3_sum: ion intensity values detected in each sample (Z1, Z2... K3)

# 3 Data Visualization 

# This code visualizes lipidomics data in the following ways
## 1. Principal component analysis of variation between samples
## 2. Ridge Plots
## 3. Scatter plot
## 4. Scree plot
## 5. QC Bar plot
## 6. Bar plot per class
## 7. Heatmaps


## R Package Versions

The versions of the R packages used in this project are as follows:

- **dplyr**: 1.1.2
- **tidyverse**: 2.0.0
- **ggplot2**: 3.4.2
- **ggrepel**: 0.9.3
- **ggridges**: 0.5.4
- **pheatmap**: 1.0.12
- **reshape**: 0.8.9
- **readxl**: 1.4.2
- **readr**: 2.1.4
- **tidyr**: 1.3.0
- **edgeR**: 3.42.4