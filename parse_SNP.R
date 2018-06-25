parse_SNP <- function(all_data, LD) {

  # give upstream and downstream their own variables to make life easier
  LD_upstream <- LD[[1]]
  LD_downstream <- LD[[2]]

  # upstream loop
  for (name in names(LD_upstream)){
    temp_data <- LD_upstream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>% 
      mutate(Marker2 = paste0("S", Locus1, "_", Position2))
    
    # filter based on SNPs in stats/effects
    temp_data <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
    
    # retrieve linked SNPs
    chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= 0.8)
    
    # catch everything that didn't make it through the filter
    chr_unlinked <- temp_data %>% filter(R.2 >= 0.8) %>% arrange(Position1, Dist_bp)
    
    # find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
    
    block_SNPs <- chr_linked %>% group_by(Marker1) %>% summarise(count = n()) %>% filter(count > 1)
    block_SNPs <- chr_linked %>% filter(Marker1 %in% block_SNPs$Marker1)
    index <- c(0, cumsum(abs(diff(block_SNPs$Site2)) > 1))
    blocks <- split(block_SNPs, index)
  
    test_data <- merge(block_SNPs, all_data, by.x = "Marker1", by.y = "Marker") %>%
      select(Locus1, Position1, Position2, Site2, Dist_bp, R.2, p, Estimate.x, Marker2)
    test_data <- merge(test_data, all_data, by.x = "Marker2", by.y = "Marker")
    
    test_data <- mutate(test_data, SNP1_pval = test_data$p)
    test_data <- mutate(test_data, SNP1_effect = test_data$Estimate.x)
    single_SNPs <- chr_linked %>% group_by(Position1) %>% summarise(count = n()) %>% filter(count == 1)
    
    ## Last two steps are below this line. All your stuff should be above this.
    
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
