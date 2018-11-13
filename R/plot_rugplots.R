plot_pathways <- function(rugplots_data, filter_type, filter_parameter, mode, output_directory) {
  
  rugplots_data <- rugplots_data %>% arrange(NESrank)
  write.table(rugplots_data, file=paste0(output_directory, "/", mode, ".full.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  
  if (filter_type == "rank") {
    rugplots_data <- rugplots_data %>% filter(NESrank <= filter_parameter)
  } else if (filter_type == "pvalue") {
    rugplots_data <- rugplots_data %>% filter(pvalue <= filter_parameter)
  } else if (filter_type == "FDR") {
    rugplots_data <- rugplots_data %>% filter(qvalue <= filter_parameter)
  } else {
    print("Incorrect filtering type. Filtering a p-value <= 0.05")
    rugplots_data <- rugplots_data %>% filter(pvalue <= 0.05)
  }
  
  write.table(rugplots_data, file=paste0(output_directory, "/", mode, ".filtered.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  rugplots_split <- split(rugplots_data, rugplots_data$NESrank)
  
  for (rank in names(rugplots_split)) {
    temp_data <- rugplots_split[[rank]]
    title = paste0(unique(as.character(temp_data$pathway_id)), " - ", unique(as.character(temp_data$pathway_name)))
    NES_Observed = unique(as.character(temp_data$NES_Observed))
    intercept <- temp_data %>% arrange(desc(phits_pmisses)) %>% select(rank)
    intercept <- intercept[,1][1]
    rugplot <- ggplot(temp_data, aes(x = rank, y = phits_pmisses)) +
      geom_line(stat = "identity") +
      geom_rug(sides = "t", position = "jitter") + 
      geom_vline(xintercept = intercept, color = "black", linetype = "longdash") + 
      ggtitle(title) + 
      labs(x = "Gene Rank", y = "Running Enrichment Score") + 
      scale_x_continuous(breaks = c(0,5000,10000,15000,20000,25000)) + 
      theme(axis.text = element_text (color = "black"), panel.background = element_rect (color = "black", fill = "pink"))
    ggsave(paste0(output_directory, "/", mode, ".", unique(as.character(temp_data$pathway_id)), ".png"), rugplot)
  }
}
  