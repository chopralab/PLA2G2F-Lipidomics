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

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(readxl)
library(readr)
library(tidyr)

# set WD
getwd()
setwd("")

# read in raw data
filepath = "../data/" # make sure raw data file is placed in this folder
raw_coltypes = c(rep("text", 6), "numeric", "numeric", "text", "text", 
                 rep("numeric",4), "text", rep("numeric", 9)) # specify column types so scientific notation gets read as numbers
raw_df <- read_xls(paste0(filepath, "20220613_Cinzia Control and cKO.xls"), col_types = raw_coltypes) %>%
    mutate(across(blank:K3, ~ replace_na(., 0))) # convert NaN to 0

# sum lipids with identical formulas
new_df <- raw_df %>% 
    group_by(LipidMolec) %>%
    transmute(Z1_sum = sum(Z1), Z2_sum = sum(Z2), Z3_sum = sum(Z3), Z4_sum = sum(Z4), Z5_sum = sum(Z5), 
              K1_sum = sum(K1), K2_sum = sum(K2), K3_sum = sum(K3), blank_sum = sum(blank)) %>%
    distinct(LipidMolec, .keep_all = TRUE) # this removes other text columns so need to add them back

# get lipid IDs from raw data
IDs <- raw_df %>% 
    select(LipidMolec:TotalGrade) %>%
    distinct(LipidMolec, .keep_all = TRUE)

# put IDs back
results_all <- merge(IDs, new_df, by = "LipidMolec")
results_AB <- filter(results_all, TotalGrade != "C")

# write results
output_filepath = "../filtered_data/"
dir.create(output_filepath, F)
write_csv(results_all, paste0(output_filepath, "sum_identical_lipids_all.csv")) # 1493 unique lipids
write_csv(results_AB, paste0(output_filepath, "sum_identical_lipids_AB_only.csv")) # 736 grade A or B lipids

# table for all detected lipids
table(results_all$ClassKey) 
table_all <- arrange(as.data.frame(table(results_all$type)), Freq)
write_csv(table_all, file = paste0(output_filepath, "/all_lipids_detected.csv"))

table(results_AB$ClassKey) 
table_AB <- arrange(as.data.frame(table(results_AB$type)), Freq)
write_csv(table_AB, file = paste0(output_filepath, "/grade_AB_lipids_detected.csv")) # this is then visualized in Graphpad Prism

# sum ion intensity for each lipid class
sum_all <- results_all %>%
    group_by(ClassKey) %>%
    summarise(Z1_TotalIonIntensity = sum(Z1_sum), Z2_TotalIonIntensity = sum(Z2_sum), Z3_TotalIonIntensity = sum(Z3_sum),
              Z4_TotalIonIntensity = sum(Z4_sum), Z5_TotalIonIntensity = sum(Z5_sum), K1_TotalIonIntensity = sum(K1_sum),
              K2_TotalIonIntensity = sum(K2_sum), K3_TotalIonIntensity = sum(K3_sum))
write_csv(sum_all, file = paste0(output_filepath, "/sum_ion_intensity_per_class_all.csv"))

sum_AB <- results_AB %>%
    group_by(ClassKey) %>%
    summarise(Z1_TotalIonIntensity = sum(Z1_sum), Z2_TotalIonIntensity = sum(Z2_sum), Z3_TotalIonIntensity = sum(Z3_sum),
              Z4_TotalIonIntensity = sum(Z4_sum), Z5_TotalIonIntensity = sum(Z5_sum), K1_TotalIonIntensity = sum(K1_sum),
              K2_TotalIonIntensity = sum(K2_sum), K3_TotalIonIntensity = sum(K3_sum))
write_csv(sum_AB, file = paste0(output_filepath, "/sum_ion_intensity_per_class_grade_AB.csv")) # this is then visualized in Graphpad Prism

