read_gff <- function(gff_file) {
  
  gff <- read.table(gff_file, comment.char = "#", col.names = 
               c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))  %>% 
    select(chr, attributes, start, end) %>% 
    mutate(tempcol = stringr::str_split(attributes, ';Name=')) %>%
    rowwise() %>%
    mutate(name = gsub("ID=", "", unlist(tempcol)[1]), tempcol = NULL, attributes = NULL)
  
  split(gff, gff$chr)
}

find_genes <- function(gff, snp_df, window) {
  
  no_linked = snp_df %>% filter(linkedSNP_count == 0) %>% mutate(Position = Position1) %>% 
    select(-c(Position1, Position2, Site1, Site2, Dist_bp))
  one_linked = snp_df %>% filter(linkedSNP_count == 1) %>% mutate(Position = ifelse(SNP1_effect != SNP2_effect,
                                                                                    ifelse(SNP1_effect < 0 & SNP2_effect < 0,
                                                                                           ifelse(SNP1_effect < SNP2_effect,
                                                                                                  Position1,
                                                                                                  Position2),
                                                                                           ifelse(SNP1_effect > 0 & SNP2_effect > 0,
                                                                                                  ifelse(SNP1_effect > SNP2_effect,
                                                                                                         Position1,
                                                                                                         Position2),
                                                                                                  ifelse(SNP1_effect * SNP2_effect < 0,
                                                                                                         ifelse(SNP1_pval != SNP2_pval,
                                                                                                                ifelse(SNP1_pval < SNP2_pval,
                                                                                                                       Position1,
                                                                                                                       Position2),
                                                                                                                "problem"),
                                                                                                         Position2))),
                                                                                    Position2)) %>% 
    select(-c(Position1, Position2, Site1, Site2, Dist_bp))
  more_linked = snp_df %>% filter(linkedSNP_count > 1) %>% mutate(Position = ifelse(Dist_bp <= window,
                                                                                    ifelse(SNP1_effect == SNP2_effect,
                                                                                           Position2,
                                                                                           ifelse(SNP1_effect < 0 & SNP2_effect < 0,
                                                                                                  ifelse(SNP1_effect < SNP2_effect,
                                                                                                         Position1,
                                                                                                         Position2),
                                                                                                  ifelse(SNP1_effect > 0 & SNP2_effect > 0,
                                                                                                         ifelse(SNP1_effect > SNP2_effect,
                                                                                                                Position1,
                                                                                                                Position2),
                                                                                                         ifelse(SNP1_effect * SNP2_effect < 0,
                                                                                                                ifelse(SNP1_pval < SNP2_pval,
                                                                                                                       Position1,
                                                                                                                       Position2),
                                                                                                                "problem")))),
                                                                                    Position2)) %>%
    select(-c(Position1, Position2, Site1, Site2, Dist_bp))
  problem_linked = snp_df %>% filter(linkedSNP_count == -1) %>% mutate(Position = ifelse(SNP1_effect != SNP2_effect,
                                                                                         ifelse(SNP1_effect < 0 & SNP2_effect < 0,
                                                                                                ifelse(SNP1_effect < SNP2_effect,
                                                                                                       Position1,
                                                                                                       Position2),
                                                                                                ifelse(SNP1_effect > 0 & SNP2_effect > 0,
                                                                                                       ifelse(SNP1_effect > SNP2_effect,
                                                                                                              Position1,
                                                                                                              Position2),
                                                                                                       ifelse(SNP1_effect * SNP2_effect < 0,
                                                                                                              "problem",
                                                                                                              "problem"))),
                                                                                         Position2)) %>% 
    select(-c(Position1, Position2, Site1, Site2, Dist_bp))
  
  snp_df = rbind(no_linked, one_linked, problem_linked, more_linked) %>% 
    filter(Position != "problem") %>% mutate(Position = as.integer(Position)) %>% 
    mutate(window_start = Position-1000, window_end = Position + 1000) %>% mutate(chr = Locus1)
  inner_join(snp_df, gff, by="chr") %>% filter((window_start >= start & window_end <= end) | (start >= window_start & end <= window_end) | 
                                                 (window_start >= start & window_start <= end) | (window_end >= start & window_end <= end)) %>% 
    select(-c(chr, start, end, window_start, window_end))
}

parse_block <- function(block) {
  block <- block %>% mutate(linkedSNP_count = nrow(block))
  
  # count the number of negative/positive effects in each block
  negative <- sum(block$SNP2_effect < 0)
  positive <- sum(block$SNP2_effect > 0)
  
  # Find SNP with largest negative or positive effect
  if (positive > negative) {
    # sort in descending order
    block <- block %>% arrange(desc(SNP2_effect))
    
    # get top row
    block <- block[1,]
    
  } else if (negative > positive) {
    # sort in ascending order
    block <- block %>% arrange(SNP2_effect)
    
    # get top row
    block <- block[1,]
    
  } else if (negative == positive) {
    if(block$SNP1_effect[[1]] > 0){
      # sort in descending order
      block <- block %>% arrange(desc(SNP2_effect))
      
      # get top row
      block <- block[1,] %>% mutate(linkedSNP_count = -1)
    } else{
      # sort in ascending order
      block <- block %>% arrange(SNP2_effect)
      
      # get top row
      block <- block[1,] %>% mutate(linkedSNP_count = -1)
    }
  }
  block
}

tag_SNPs <- function(untagged_genes) {
  neg_genes <- untagged_genes %>% filter(SNP2_effect < 0) %>% arrange(SNP2_effect, SNP2_pval)
  pos_genes <- untagged_genes %>% filter(SNP2_effect > 0) %>% arrange(desc(SNP2_effect, SNP2_pval))
  negative <- sum(neg_genes$linkedSNP_count)
  positive <- sum(pos_genes$linkedSNP_count)
  
  if (positive > negative){
    untagged_genes <- untagged_genes %>% arrange(desc(SNP2_effect, SNP2_pval))
    tagSNP <- untagged_genes[1,]
  } else if(negative > positive){
    untagged_genes <- untagged_genes %>% arrange(SNP2_effect, SNP2_pval)
    tagSNP <- untagged_genes[1,]
  } else if(positive == negative){
    pos_max <- pos_genes[1,]
    neg_max <- neg_genes[1,]
    neg_max$SNP2_effect <- abs(neg_max$SNP2_effect)
    
    if(pos_max$SNP2_effect > neg_max$SNP2_effect){
      tagSNP <- pos_max
    }
    else{
      tagSNP <- neg_max
    }
  }
  tagSNP %>% mutate(linkedSNP_count=negative+positive)
}

parse_SNP <- function(all_data, LD, gff_file, window, r_squared_cutoff, num_cores) {
  
  full_gff <- read_gff(gff_file)
  
  cl <- parallel::makeCluster(num_cores, outfile="")
  registerDoParallel()
  
  # UP/DOWNSTREAM LOOP
  for (i in 1:length(LD)) {
    LD_stream <- LD[[i]]

    # BEGIN PROCESSING BY CHROMOSOMES LOOP
    for (name in names(LD_stream)) {
      temp_data <- LD_stream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>%
        mutate(Marker2 = paste0("S", Locus1, "_", Position2))
      
      chr_data <- all_data %>% filter(Chr == as.integer(name))
  
      # retrieve linked SNPs
      chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= r_squared_cutoff)
      
      # get block SNPs
      block_SNPs <- chr_linked %>% group_by(Marker1) %>% dplyr::summarise(count = n()) %>% filter(count > 1)
      block_SNPs <- chr_linked %>% filter(Marker1 %in% block_SNPs$Marker1)
      
      # look up p-value and effect data for SNP1
      block_SNPs <- merge(block_SNPs, chr_data, by.x = "Marker1", by.y = "Marker") %>% 
        mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2)
      
      # look up p-value and effect data for SNP2
      block_SNPs <- merge(block_SNPs, chr_data, by.x = "Marker2", by.y = "Marker") %>%
        mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect) %>% 
        arrange(Position1)
      
      # split into blocks by Position1 and difference in Site2 where difference > 1
      index <- c(0, cumsum(abs(diff(block_SNPs$Site2)) > 1))
      block_SNPs_list <- split(block_SNPs, paste(block_SNPs$Position1, index))
      
      block_SNPs <- foreach(block_SNPs_item=block_SNPs_list, .combine = rbind, .packages = c('dplyr', 'past')) %dopar% {
          parse_block(block_SNPs_item)
      }
      
      # get all the single SNPs and get effects/p-value data
      single_SNPs <- chr_linked %>% group_by(Marker1) %>% dplyr::summarise(count = n()) %>% filter(count == 1)
      single_SNPs <- chr_linked %>% filter(Marker1 %in% single_SNPs$Marker1) %>% mutate(linkedSNP_count = 1)
      single_SNPs <- merge(single_SNPs, chr_data, by.x = "Marker1", by.y = "Marker") %>% 
        mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)
      single_SNPs <- merge(single_SNPs, chr_data, by.x = "Marker2", by.y = "Marker") %>%
        mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count)
      
      # get all unlinked SNP data
      chr_unlinked <- temp_data %>% filter(R.2 < r_squared_cutoff) %>% arrange(Position1, Dist_bp)
      
      # find the first instance of every unlinked SNP
      unlinked_SNPs <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]
      
      # get all the unlinked SNPs and get effects/p-value data
      unlinked_SNPs <- unlinked_SNPs %>% mutate(linkedSNP_count = 0)
      unlinked_SNPs <- merge(unlinked_SNPs, chr_data, by.x = "Marker1", by.y = "Marker") %>% 
        mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)
      unlinked_SNPs <- merge(unlinked_SNPs, chr_data, by.x = "Marker2", by.y = "Marker") %>%
        mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count)
      
      # combined block tagSNPs, single_SNPs, unlinked SNPs
      # each line is independent so this data can be split
      # for processing in parallel
      SNPs = rbind(block_SNPs, single_SNPs, unlinked_SNPs)
      
      # get positions and split them
      split = num_cores*2
      SNPS_list <- split(SNPs, rep(1:split, length.out = nrow(SNPs), each = ceiling(nrow(SNPs)/split)))
      
      # initialize genes dataframe
      SNP_genes = NULL
      
      # subset gff to only handle this chromosome
      gff <- full_gff[[name]] %>% mutate(chr = as.integer(chr))
      
      # get genes in parallel
      all_genes <- foreach(SNPs_item=SNPS_list, .combine = rbind, .packages = c('dplyr', 'past')) %dopar%{
          find_genes(gff, SNPs_item, window)
      }
      
      # store modified and filtered data
      LD_stream[[name]] <- all_genes
      
    } # END UP/DOWNSTREAM LOOP
    
    # set up/downstream to modified data
    LD[[i]] <- LD_stream
    
  }
  
  combined_streams = NULL
  for (name in names(LD[[1]])) {
    all_genes <- rbind(LD[[1]][[name]], LD[[2]][[name]]) %>% 
      arrange(name) %>% 
      mutate(linkedSNP_count = ifelse(linkedSNP_count <= 0, 1 ,linkedSNP_count))
    group_genes <- split(all_genes, f = all_genes$name)
    tagged_genes <- foreach(block=group_genes, .combine = rbind, .packages = c('dplyr', 'past')) %dopar%{
      tag_SNPs(block)
    }
    combined_streams[[name]] <- tagged_genes
  }
  
  stopCluster(cl)
  
  # combined genes from each chromosome into single dataframe
  genes = NULL
  for (name in names(combined_streams)) {
    genes <- rbind(genes, combined_streams[[name]])
  }
  
  # return list of genes
  genes
}
