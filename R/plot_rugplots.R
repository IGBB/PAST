plot_pathways <- function(rugplots_data, output_directory) {
  
  rugplots_data <- rugplots_data %>% arrange(NESrank)
  rugplots_split <- split(rugplots_data, rugplots_data$NESrank)
  
  for (rank in names(rugplots_split)) {
    temp_data <- rugplots_split[[rank]]
    title = paste0(unique(as.character(temp_data$pathway_id)), " - ", unique(as.character(temp_data$pathway_name)))
    NES_Observed = unique(as.character(temp_data$NES_Observed))
    intercept <- temp_data %>% arrange(desc(phit_pmiss)) %>% select(rank)
    intercept <- intercept[,1][1]
    rugplot <- ggplot(temp_data, aes(x = rank, y = phit_pmiss)) +
      geom_line(stat = "identity") +
      geom_rug(sides = "t", position = "jitter") + 
      geom_vline(xintercept = intercept, color = "black", linetype = "longdash") + 
      ggtitle(title) + 
      labs(x = "Gene Rank", y = "Running Enrichment Score") + 
      scale_x_continuous(breaks = c(0,5000,10000,15000,20000,25000)) + 
      theme(axis.text = element_text (color = "black"), panel.background = element_rect (color = "black", fill = "pink"))
    ggsave(paste0(output_directory, "/", unique(as.character(temp_data$pathway_id)), ".png"), rugplot)
  }
}
  