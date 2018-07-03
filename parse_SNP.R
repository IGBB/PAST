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
    
    # find all unique positions in block_SNPs
    positions <- block_SNPs %>% group_by(Position1) %>% summarise(count = n())
   
    
   
     
     # begin block SNP loops
    for(pos in positions$Position1){
       positions_block_SNPs <- block_SNPs %>% filter(Position1 == pos)
       index <- c(0, cumsum(abs(diff(positions_block_SNPs$Site2)) > 1))
       blocks <- split(positions_block_SNPs, index)
      for (block_name in names(blocks)){
        block <- blocks[[block_name]]
        
        #I moved the merge command back to this spot <<<<<<<<
        # merge block and positions by Position1 data and rename count to linkedSNP_count
        block <- block %>% merge(block_SNPs, positions, by = "Position1") %>%
        mutate(linkedSNP_count = nrow(block))
      
        # merge the block SNPs with the data in all data to find the data for the markers in block SNPs
        block <- merge(block, all_data, by.x = "Marker1", by.y = "Marker") %>% 
          mutate(SNP1_pval = p, SNP1_effect = Estimate.x) %>% 
          select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)
      
        # merge the block SNPs with all data to get the marker data for Marker2 SNPs
        block <- merge(block, all_data, by.x = "Marker2", by.y = "Marker") %>%
          mutate(SNP2_pval = p, SNP2_effect = Estimate.x) %>%
          select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count)
        
        # count the number of negative/positive effects in each week
        temp_stuff <- block[[block_name]]
        negative <- sum(temp_stuff$SNP1_effect < 0)
        positive <- sum(temp_stuff$SNP1_effect > 0)
        
        if (positive > negative){ #Positive
          print("Positive")
          
          #sort in descending order
          order(-temp_stuff$SNP1_effect)
          temp_stuff <- block[[block_name[[1]]]]
        } else if (negative > positive){ #Negative
          print("Negative")
          
          # sort in ascending order
          order(temp_stuff$SNP1_effect)
        } else if (negative == positive){ #Equal
          print ("Equal")
        }
      }
    }
  }
     
    
   # get all the single SNPs in a data frame by themselves
    single_SNPs <- chr_linked %>% group_by(Position1) %>% summarise(count = n()) %>% filter(count == 1)
    
    ## Last two steps are below this line. All your stuff should be above this.
    
    # combine linked and unlinked and arrange Position1 in ascending order
    temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
    
    # store modified and filtered data
    LD_upstream[[name]] <- temp_data
  }

  # downstream loop
  for (name in names(LD_downstream)){
    temp_data <- LD_downstream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>%
      mutate(Marker2 = paste0("S", Locus1, "_", Position2))
    
    # filter based on SNPs in stats/effects
    temp_data <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
    
    # retrieve linked SNPs
    chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= 0.8)
    
    # catch everything that didn't make it through the filter
    chr_unlinked <- temp_data %>% filter(R.2 >= 0.8) %>% arrange(Position1)
    
    # find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
    
    
    block_SNPs <- chr_linked %>% group_by(Marker1) %>% summarise(count = n()) %>% filter(count > 1)
    block_SNPs <- chr_linked %>% filter(Marker1 %in% block_SNPs$Marker1)
    
    # merge the block SNPs with the data in all data to find the data for the markers in block SNPs
    block_SNPs <- merge(block_SNPs, all_data, by.x = "Marker1", by.y = "Marker") %>% 
      mutate(SNP2_pval = p, SNP2_effect = Estimate.x) %>% 
      select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP2_pval, SNP2_effect, Marker2)
    
    # merge the block SNPs with all data to get the marker data for Marker2 SNPs
    block_SNPs <- merge(block_SNPs, all_data, by.x = "Marker2", by.y = "Marker") %>%
      mutate(SNP1_pval = p, SNP1_effect = Estimate.x) %>%
      select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect)
    
    index <- c(0, cumsum(abs(diff(block_SNPs$Site2)) > 1))
    blocks <- split(block_SNPs, index)
   
    for (stuff in names(blocks)){
      temp_stuff <- block[[stuff]]
      negative <- sum(temp_stuff$SNP1_effect < 0)
      positive <- sum(temp_stuff$SNP1_effect > 0)
      if (positive > negative){ #Positive
        print("Positive")
        
        # put in descending order
        order(-temp_stuff$SNP1_effect)
      } else if (negative > positive){ #Negative
        print("Negative")
        
        # put in ascending order
        order(temp_stuff$SNP1_effect)
        temp_stuff <- block[[stuff[[1]]]]
      } else if (negative == positive){ #Equal
        print ("Equal")
      }
      
      
    }
    
    # get all the single SNPs in a data frame by themselves
    single_SNPs <- chr_linked %>% group_by(Position1) %>% summarise(count = n()) %>% filter(count == 1)
    
    # combine linked and unlinked and arrange Position1 in ascending order
    temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
    
    # store modified and filtered data
    LD_downstream[[name]] <- temp_data
  }
  
  # return modified lists
  list(LD_upstream, LD_downstream)
}
