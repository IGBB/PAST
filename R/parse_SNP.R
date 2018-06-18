parse_SNP <- function(all_data, LD) {

  #Give upstream and downstream their own variables to make life easier
  LD_upstream <- LD[[1]]
  LD_downstream <- LD[[2]]

  for (name in names(LD_upstream)){
    temp_data <- LD_upstream[[name]]
    LD_upstream[[name]] <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
  }

  for (name in names(LD_downstream)){
    temp_data <- LD_downstream[[name]]
    LD_downstream[[name]] <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
  }
  
  names(LD_upstream)
  for (name in names(LD_upstream))
  {
    print(name)
    print(LD_upstream[[1]]$Position1)
  }
  
  for (chr in names(LD_upstream)){
    chr_up <- LD_upstream[[chr]]
    
    #Retrieve linked SNPs
    chr_linked <- chr_up %>% arrange(Position1) %>% filter(R.2 >= 0.8)
    
    #Catch everything that didn't make it through the filter
    chr_unlinked <- chr_up %>% filter(R.2 >= 0.8) %>% arrange(Position1, Dist_bp)
    
    #Find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
    
    #Combine linked and unlinked and arrange Position in ascending order
    chr_up <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
  }
  
  #Same steps as the loop prior but for LD_downstream
  for (chr in names(LD_downstream)){
    chr_down <- LD_downstream[[chr]]
    chr_linked <- chr_down %>% arrange(Position1) %>% filter(R.2 >= 0.8)
    
    #Only arrange by Position1 not Position1 and Dist_bp
    chr_unlinked <- chr_down %>% filter(R.2 >= 0.8) %>% arrange(Position1)
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
    chr_down <-  rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
  }
}
