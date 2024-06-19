# This code performs EdgeR analysis for pre-processed data.
## EdgeR is a form of differential expression analysis which uses a generalized linear model (GLM). 
## The GLM uses a negative binomial distribution to model lipid ion count data (or mRNA count data) accounting for the technical and biological variability. 

# References
citation("edgeR")

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

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(readxl)
library(edgeR)
library(tidyverse)

# set WD
getwd()
setwd("")

# read in pre-processed data
file_path = "../filtered_data/" # make sure filtered data are in this folder
file_list = list.files(path=file_path, pattern=NULL, all.files=FALSE, full.names=FALSE)
file_list

read_lipids <- function(file_path, filename) {
    data<- read_csv(paste0(file_path, filename), show_col_types = FALSE) %>%
    dplyr::rename(lipid = 'LipidMolec', type= 'ClassKey') %>%
    mutate_if(is.numeric, function(x) replace(x, is.nan(x), 0)) # replace NaN values with 0
}

data_all <- read_lipids(file_path, "sum_identical_lipids_all.csv")
data_AB <- read_lipids(file_path, "sum_identical_lipids_AB_only.csv")

# set EdgeR groups
## for this experiment, comparing Pla2g2f KO vs WT
gr1_name = "WT" # Z group = WT
gr1_length = 5 # number of samples in Z group
gr1 = rep(gr1_name, gr1_length)

gr2_name = "KO" # K group = Pla2g2f KO
gr2_length = 3 # number of samples in K group
gr2 = rep(gr2_name, gr2_length)

blank_name = "blank" # column name for blank

groups_expr1 = c(gr1, gr2, blank_name) %>% # create factor from groups
              factor(levels = c(blank_name, gr1_name, gr2_name)) 
design_expr1 = model.matrix(~groups_expr1) # create design matrix with experimental groups as column names
colnames(design_expr1) = c("Intercept", gr1_name, gr2_name) 
contrast_expr1 = makeContrasts( # create contrast to compare Pla2g2f KO vs WT
    H = KO - WT, 
    levels = groups_expr1
)

# EdgeR functions
perform_raw_analysis <- function(data, counts_col, lipid_col, group, design) {
    y <- DGEList(counts = data[,counts_col],
                 genes = data[,lipid_col], 
                 group = group) # create DGEList object from data
    y <- calcNormFactors(y, method = "TMM") # calculate scaling factor for raw data
    y <- estimateCommonDisp(y, design = design) # estimate a common dispersion value across all genes/lipids
    y
}

calculate_significance <- function(data, design, contrast){
    data %>% 
    glmFit(design = design) %>% # fit NB model to the counts for each gene/lipid
    glmLRT(contrast = contrast) # conduct likelihood ratio tests for the given contrasts (to be tested equal to 0)
}

# perform EdgeR analysis
colnames(data_AB)
counts_cols = 16:24 # where ion intensities are in the raw data frame
lipid_cols = c(1,3) # where lipid names & types are in the raw data frame

results_all <- data_all %>%
    perform_raw_analysis(counts_cols, lipid_cols, groups_expr1, design_expr1) %>% 
    calculate_significance(design_expr1, contrast_expr1) %>% 
    topTags(5000) %>% 
    as.data.frame() %>% 
    remove_rownames

results_AB <- data_AB %>%
    perform_raw_analysis(counts_cols, lipid_cols, groups_expr1, design_expr1) %>% 
    calculate_significance(design_expr1, contrast_expr1) %>% 
    topTags(5000)  %>% 
    as.data.frame() %>% 
    remove_rownames

# combine original lipid intensities & EdgeR outputs
results_all <- merge(data_all[,c(1, counts_cols)], results_all, by = 'lipid')
results_AB <- merge(data_AB[,c(1, counts_cols)], results_AB, by = 'lipid')

# add group mean to results dataframe
add_means <- function(data, g1_cols, g2_cols){ # this function calculate means for 2 groups and add them to the dataframe
    data <- data %>%
    mutate(mean1 = rowMeans(select(., all_of(g1_cols))),
           mean2 = rowMeans(select(., all_of(g2_cols)))) # calculate row means for each group
}

head(results_AB)
grp1_cols <- 2:6 # specify where group 1 & 2 columns are in the dataframe
grp2_cols <- 7:9

data_all_means <- add_means(results_AB, grp1_cols, grp2_cols) %>%
    select(lipid, type, mean1, mean2, logFC, logCPM, LR, PValue, FDR, Z1_sum:K3_sum)
data_AB_means <- add_means(results_AB, grp1_cols, grp2_cols) %>%
    select(lipid, type, mean1, mean2, logFC, logCPM, LR, PValue, FDR, Z1_sum:K3_sum)

# export results
out_filepath = "../edger_results/"
dir.create(path = out_filepath, F)
write_csv(data_all_means, file = paste0(out_filepath, 'Pla2g2fko_vs_WT_all.csv'))
write_csv(data_AB_means, file = paste0(out_filepath, 'Pla2g2fko_vs_WT_gradeAB.csv'))

