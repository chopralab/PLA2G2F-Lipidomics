# Code for functions used to create PCA plots

# Calculate principal components
calculate_pc <- function(data, counts_col) {
    numerical_data <- t(data[,counts_col]) # transpose data so that row contain observations/groups, column contains variables
    numerical_data <- numerical_data[, apply(numerical_data, 2, function(x) any(x != x[1]))] # remove columns with identical values
    prcomp(numerical_data, center = TRUE, scale = TRUE) # returns PCA results
}

# Scree plot
make_scree_plot <- function(pca, title) {
    variance = pca$sdev^2 / sum(pca$sdev^2) # calculate total variance experienced by each pc
    as.data.frame(variance) %>%  # convert vector to tibble for plotting 
        mutate(PC = 1:length(variance)) %>%
        ggplot() + 
        geom_col(aes(x = PC, y = variance))+
        xlab("Principal Component") + 
        ylab("Variance") +
        ggtitle(title) +
        ylim(0, 1) + 
        theme_classic(base_size = 20) +
        theme(aspect.ratio = 1)
}

# PCA plot
make_pca_plot <- function(pca_result, groups, title) {
    pca_scores <- as.data.frame(pca_result$x) # get dataframe from pca results
    plot_data <- data.frame(PC1 = pca_scores$PC1,
                            PC2 = pca_scores$PC2,
                            group = groups, 
                            sample = str_sub(rownames(pca_scores), 1, 2))
    plot_data %>%
        ggplot(aes(x = PC1, y = PC2, color = group, fill = group)) +
        geom_point(size = 3) +
        geom_text_repel(aes(label = sample), size = 10) +
        stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.2)  +
        theme_classic(base_size = 20) +
        labs(color = "group") +
        ggtitle(paste0("PCA for ", title)) +
        scale_color_manual(values = c("red", "black")) +
        scale_fill_manual(values = c("red", "black")) +
        xlab(paste0("PC1: ", round((pca_result$sdev[1]^2 / sum(pca_result$sdev^2)) * 100, 2), "% variance")) +
        ylab(paste0("PC2: ", round((pca_result$sdev[2]^2 / sum(pca_result$sdev^2)) * 100, 2), "% variance")) +
        theme(aspect.ratio = 1)
}