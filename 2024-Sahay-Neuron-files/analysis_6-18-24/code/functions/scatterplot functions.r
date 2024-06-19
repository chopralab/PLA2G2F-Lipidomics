# Function to make scatterplots

make_scatter <- function(data, title, n){
    data_sig <- filter(data, FDR<0.1) %>% arrange(FDR)
    plot <- data %>%
        ggplot(aes(x = log10(mean1), y = log10(mean2))) +
        geom_point(color = "grey") + # plot all points in grey
        geom_point(data = data_sig, aes(color = logFC)) + # plot FDR<0.1 points in LogFC color
        geom_text_repel(data=head(data_sig, n),
                        aes(label=lipid),  
                        size=5,
                        nudge_y = 0.2,    
                        nudge_x = 0.2) +
        xlab("Log10(WT/Blank)") +
        ylab("Log10(KO/Blank)") +
        ggtitle(paste0("Scatter plot for ", title)) +
        theme_classic(base_size = 20) +
        theme(aspect.ratio = 1) +
        scale_color_gradient(low = "blue", high = "red", guide = "colourbar")
    plot
}
