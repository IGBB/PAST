parse_SNP <- function(all_data, LD) {

  # give upstream and downstream their own variables to make life easier
  LD_upstream <- LD[[1]]
  LD_downstream <- LD[[2]]

  # upstream loop
  for (name in names(LD_upstream)){
    temp_data <- LD_upstream[[name]]
    
    # filter based on SNPs in stats/effects
    temp_data <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
    
    # retrieve linked SNPs
    chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= 0.8)
    
    # catch everything that didn't make it through the filter
    chr_unlinked <- temp_data %>% filter(R.2 >= 0.8) %>% arrange(Position1, Dist_bp)
    
    # find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
    
    # combine linked and unlinked and arrange Position1 in ascending order
    temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
    
    # store modified and filtered data
    LD_upstream[[name]] <- temp_data
  }

  # downstream loop
  for (name in names(LD_downstream)){
    temp_data <- LD_downstream[[name]]
    
    # filter based on SNPs in stats/effects
    temp_data <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
    
    # retrieve linked SNPs
    chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= 0.8)
    
    # catch everything that didn't make it through the filter
    chr_unlinked <- temp_data %>% filter(R.2 >= 0.8) %>% arrange(Position1)
    
    # find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
    
    # combine linked and unlinked and arrange Position1 in ascending order
    temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
    
    # store modified and filtered data
    LD_downstream[[name]] <- temp_data
  }
  
  # return modified lists
  list(LD_upstream, LD_downstream)
}
