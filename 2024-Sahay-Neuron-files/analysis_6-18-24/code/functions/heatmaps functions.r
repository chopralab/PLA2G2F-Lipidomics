# Heatmaps

## calculate z scores for heatmaps
calculate_z_scores <- function(data, FDR_cutoff=0, cols) {
    if (FDR_cutoff==0) {
       result <- data %>%
       column_to_rownames('lipid') %>%
       arrange(type, logFC) %>%
       select(cols-1) %>% 
       apply(1, scale) %>% # scale function calculates z-scores using (x - mean(x)) / sd(x)
       t()
    } else {
       result <- data %>%
       filter(FDR<FDR_cutoff) %>% # filter for FDR < cutoff value
       column_to_rownames('lipid') %>%
       arrange(type, logFC) %>%
       select(cols-1) %>% 
       apply(1, scale) %>% 
       t()
    }
    colnames(result) <- colnames(data)[2:9]
    result
}

## make heatmaps
make_heatmap <- function(data, title, lipid_annotation){
    pheatmap(data, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE, 
    show_rownames = FALSE,
    border_color = NA, 
    fontsize = 10,
    legend_breaks = -2:2,
    annotation_row = lipid_annotation, 
    annotation_col = sample_ann, 
    annotation_colors = ann_colors,
    main = paste0("Heatmap for ", title), 
    filename = paste0(heatmaps_path, "Heatmap for ", title, ".pdf"))
}

make_heatmap2 <- function(data, title, lipid_annotation){
    pheatmap(data, 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE, 
    show_rownames = TRUE, # show lipid names here
    border_color = NA, 
    fontsize = 10,
    legend_breaks = -2:2,
    annotation_row = lipid_annotation, 
    annotation_col = sample_ann, 
    annotation_colors = ann_colors,
    main = paste0("Heatmap for ", title), 
    filename = paste0(heatmaps_path, "Heatmap for ", title, ".pdf"))
}
