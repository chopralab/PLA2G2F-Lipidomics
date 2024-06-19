# Functions to make ridge plot

ridge_plot <- function(data, title) {
    data %>%
    ggplot(aes(x = logFC, y = type, fill = type)) +
    geom_density_ridges2(alpha = 0.5, size = 0) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggtitle(paste0("Ridge plot for ", title)) +
    xlab("Log2(Fold Change)") +
    ylab("") +
    scale_alpha(guide = 'none') +
    scale_fill_manual(values = lipid_colors, guide = 'none') +
    theme(aspect.ratio = 2) +
    theme_classic(base_size = 20)
}