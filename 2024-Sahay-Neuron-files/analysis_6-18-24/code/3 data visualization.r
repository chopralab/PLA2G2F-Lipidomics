# This code visualizes lipidomics data in the following ways
## 1. Principal component analysis of variation between samples
## 2. 

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(pheatmap)
library(reshape)

# # Function to get package version
# get_package_version <- function(package_name) {
#   package_version <- as.character(packageVersion(package_name))
#   return(package_version)
# }
# 
# # List of packages
# packages <- c("dplyr", "tidyverse", "ggplot2", "ggrepel", "ggridges", "pheatmap", "reshape")
# 
# # Get versions
# package_versions <- sapply(packages, get_package_version)
# 
# # Print package versions
# print(package_versions)

# set WD
getwd()
setwd("")

# read in pre-processed & analyzed data
processed_filepath = "../filtered_data/"
processed_file_list = list.files(path=processed_filepath, pattern=NULL, all.files=FALSE, full.names=FALSE)
processed_file_list
sum_AB <- read_csv(paste0(processed_filepath, "sum_ion_intensity_per_class_grade_AB.csv"), show_col_types = FALSE)

analyzed_filepath = "../edgeR_results/"
analyzed_file_list = list.files(path=analyzed_filepath, pattern=NULL, all.files=FALSE, full.names=FALSE)
analyzed_file_list
data_AB <- read_csv(paste0(analyzed_filepath, "Pla2g2fko_vs_WT_gradeAB.csv"), show_col_types = FALSE)

# output folder
out_filepath = "../plots/"
dir.create(out_filepath, F)

# functions folder
functions_filepath = "functions/"

# title prefix for plots
title <- "736 grade A&B lipids"

# lipid colors
## color palette was generated using https://medialab.github.io/iwanthue/
lipid_classes <- unique(data_AB$type)
length(lipid_classes) # need 19 colors
lipid_colors <- c("Cer"="#9e4a6a","Ch-D7"="#55b74d","Ch"="#bf4fb4","ChE"="#a1b534","DG"="#7c63d1",
                  "Hex1Cer"="#d89c33","Hex2Cer"="#6b7ac3","LPC"="#d2612e",
                  "LPE"="#4facd8","LPI"="#cc4045","MG"="#4ebd9d","PA"="#d94a81","PC"="#397e4a",
                  "PE"="#d18ac7","PG"="#88b56d","PI"="#dc8572","PS"="#677729","SM"="#94642d","TG"="#bba65e")

# bar plots for QC
## These bar plots show total ion intensity from each class in a stack for each sampmle
## If sample was preweighed/measured before analysis, expect similar total ion intensity across samples
## If samples are similar in nature (e.g. same tissue type), expect similar distribution of lipid classes
colnames(sum_AB) <- c("type", "Z1", "Z2", "Z3", "Z4", "Z5", "K1", "K2", "K3")
sum_AB %>% 
    as.data.frame() %>%
    melt(id.vars = "type") %>% # reshape data into the long format
    ggplot(aes(x = variable, y = value, fill = type)) +
    geom_col() +
    theme_classic(base_size = 20) +
    labs(color = "group") +
    ggtitle(paste0("Breakdown of lipid classes across samples (", title, ")")) +
    xlab("Sample") +
    ylab("Total Ion Intensity") +
    theme(aspect.ratio = 1) +
    scale_fill_manual(values = lipid_colors)
# Save the plot with specified size
ggsave(filename = paste0(out_filepath, "QC bar plot for ", title, ".pdf"), width = 10, height = 10)

# PCA
source(paste0(functions_filepath, "PCA plot functions.r")) # load functions for PCA analysis & visualization

## calculate principal component
head(data_AB)
counts_col = 10:17 # where ion intensity values are in the dataframe
pca_results_AB <- calculate_pc(data_AB, counts_col)

## make scree plot
make_scree_plot(pca_results_AB, paste0("Scree plot for ", title))
ggsave(filename = paste0(out_filepath, "Scree plot for ", title, ".pdf"), width = 10, height = 10)

filter(data_AB, FDR<0.1) %>% # Scree plot with FDR<0.1 cut off
    calculate_pc(counts_col) %>%
    make_scree_plot(paste0("Scree plot for ", title, " (FDR<0.1)"))
ggsave(filename = paste0(out_filepath, "Scree plot for ", title, " FDR_0.1.pdf"), width = 10, height = 10)

## make PCA plot
grp1 = "WT"
grp2 = "Pla2g2f KO"
groups = c(rep(grp1, 5), rep(grp2, 3)) # create groups for PCA visualization

make_pca_plot(pca_results_AB, groups, title)
ggsave(filename = paste0(out_filepath, "PCA for ", title, ".pdf"), width = 10, height = 10)

filter(data_AB, FDR<0.1) %>% # PCA plot for lipids with FDR<0.1 cut off
    calculate_pc(counts_col) %>%
    make_pca_plot(groups, paste0(title, " (FDR < 0.1)"))
ggsave(filename = paste0(out_filepath, "PCA for ", title, " FDR_0.1", ".pdf"), width = 10, height = 10)

# Ridgeplots  
## Ridgeplot shows fold change as a density distribution of lipids for each class
## x-axis is Log2(Fold change KO vs WT)
## y-axis is number of lipids (not to scale as different lipid classes has different total number of lipids)
## Density distribution cannot be generated for lipid classes with <3 lipids
source("functions/ridgeplot functions.r")

ridge_plot(data_AB, title)
ggsave(paste0(out_filepath, "Ridge plot for ", title, ".pdf"), width = 10, height = 10)

filter(data_AB, FDR<0.1) %>% ridge_plot(title)
ggsave(paste0(out_filepath, "Ridge plot for ", title, " FDR_0.1.pdf"), width = 10, height = 10)

filter(data_AB, PValue<0.05) %>% ridge_plot(title)
ggsave(paste0(out_filepath, "Ridge plot for ", title, " P_0.05.pdf"), width = 10, height = 10)

filter(data_AB, PValue<0.01) %>% ridge_plot(title)
ggsave(paste0(out_filepath, "Ridge plot for ", title, " P_0.01.pdf"), width = 10, height = 10)

# Heatmaps
source("functions/heatmaps functions.r")

## calculate z scores for heatmaps
heatmap_data_AB <- calculate_z_scores(data_AB, cols = counts_col)
heatmap_data_AB_fdr <- calculate_z_scores(data_AB, FDR_cutoff = 0.1, cols = counts_col)

## heatmap annotations
sample_ann <- data.frame(Group = rep(c("WT", "Pla2g2f KO"), c(5,3)),
                         row.names = colnames(heatmap_data_AB))
lipid_ann <- data.frame(Class = data_AB$type,
                           row.names = data_AB$lipid)
ann_colors <- list(LipidClass = lipid_colors, 
                   Group = c("WT" = "black", "Pla2g2f KO" = "red"))

## make heatmaps
heatmaps_path <- paste0(out_filepath, "heatmaps/")
dir.create(heatmaps_path, F)
make_heatmap(heatmap_data_AB, title, lipid_ann)
make_heatmap(heatmap_data_AB_fdr, paste0(title, " FDR_0.1"), lipid_ann)

## Grouping lipids according to LipidMaps.org
## Sphingolipids: Cer, Hex1Cer, Hex2Cer, SM
## Cholesterol: Ch-D7, Ch, ChE
## Glycerophospholipids: LPC, LPE, LPI, PA, PC, PE, PG, PI, PS
## Glycerolipids: TG, DG, MG

sph <- c("Cer", "Hex1Cer", "Hex2Cer", "SM")
cholesterols <- c("Ch-D7", "Ch", "ChE")
gpl <- c("PA", "PC", "PE", "PG", "PI", "PS", "LPC", "LPE", "LPI")
glyceroll <- c("TG", "DG", "MG")

## heatmap colored by groups of lipids
data_AB_wclass <- data_AB %>%
    as_tibble() %>%
    mutate(group = case_when(type %in% sph ~ 'Sphingolipid',
                             type %in% cholesterols ~ 'Cholesterol',
                             type %in% gpl ~ 'Glycerophospholipid',
                             type %in% glyceroll ~ 'Glycerolipid'))
lipid_ann2 <- data.frame(Class = data_AB_wclass$group,
                        row.names = data_AB_wclass$lipid)
ann_colors <- list(LipidClass = c("Sphingolipid"="#7aa556", "Cholesterol"="#c65999",
                                  "Glycerophospholipid"="#c96d44", "Glycerolipid"="#777acd"), 
    Group = c("WT" = "black", "Pla2g2f KO" = "red"))
heatmap_data_AB_wclass <- calculate_z_scores(data_AB_wclass, cols = counts_col)
make_heatmap(heatmap_data_AB_wclass, paste0(title, " with class"), lipid_ann2)

## separate heatmaps for each group
ann_colors <- list(LipidClass = lipid_colors, 
                   Group = c("WT" = "black", "Pla2g2f KO" = "red"))

filter(data_AB_wclass, group == "Sphingolipid") %>%
    calculate_z_scores(cols = counts_col) %>%
    make_heatmap2(paste0(title, " Sphingolipids"), lipid_ann)
filter(data_AB_wclass, group == "Cholesterol") %>%
    calculate_z_scores(cols = counts_col) %>%
    make_heatmap2(paste0(title, " Cholesterol"), lipid_ann)
filter(data_AB_wclass, group == "Glycerophospholipid") %>%
    calculate_z_scores(cols = counts_col) %>%
    make_heatmap2(paste0(title, " Glycerophospholipid"), lipid_ann)
filter(data_AB_wclass, group == "Glycerolipid") %>%
    calculate_z_scores(cols = counts_col) %>%
    make_heatmap2(paste0(title, " Glycerolipid"), lipid_ann)

dev.off() # so we can keep previewing the plots

# scatter plot
## scatterplot shows the average intensity of lipids detected in WT vs. Pla2g2fKO samples
## x-axis is Log10 transformed fold-change of WT/Blank
## y-axis is Log10 transformed fold-change of Pla2g2fKO/Blank
## colored dots are lipids with FDR < 0.1, and they are colored based on fold-change of Pla2g2fKO/WT
source("functions/scatterplot functions.r")

make_scatter(data_AB, title, 14)
ggsave(filename = paste0(out_filepath, "Scatter plot for ", title, ".pdf"), width = 10, height = 10)

# bar plots for each class
## bar plots show normalized average intensity of each lipid
## After average is calculated for each group, lipids are normalized by dividing over the larger group
barplot_outpath <- paste0(out_filepath, "/bar plots per class")
dir.create(barplot_outpath, F)

bar_data <- data_AB[c("mean1","mean2")] %>% 
    apply(1, function(x)(x/max(x))) %>% # range normalization by dividing over the maximum
    t() %>%
    as_tibble() %>%
    mutate(lipid = data_AB$lipid,
            type = data_AB$type) %>%
    as.data.frame() %>%
    melt(id = c("lipid", "type")) 

for (lipid_class in lipid_classes) {
   bar_data %>%
        filter(type == lipid_class) %>%
        ggplot(aes(x = lipid, y = value, fill = variable)) +
        geom_col(position = "dodge") +
        xlab("") +
        ylab("Normalized Intensity") +
        ggtitle(paste0("Normalized intensity for ", lipid_class)) +
        theme_classic(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_fill_manual(values = c("#9A9695", "#329666"))
    # readline(prompt="Press [enter] to proceed")
    ggsave(filename = paste0(barplot_outpath, "/normalized intensity for ", lipid_class, " in ", title, ".pdf"), width = 10, height = 10)
}
# some lipid classes have too many lipids, so their names get squished on the x-axis
# these can be fixed by changing the font size of the text in a figure editor (e.g. Affinity Designer)

