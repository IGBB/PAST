read_gff <- function(gff_file) {
  
  gff <- read.table(gff_file, comment.char = "#", col.names = 
                      c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))  %>% 
    select(chr, attributes, start, end) %>% 
    mutate(tempcol = stringr::str_split(attributes, ';Name=')) %>%
    rowwise() %>%
    mutate(name = gsub("ID=", "", unlist(tempcol)[1]), tempcol = NULL, attributes = NULL, chr=gsub("chr", "", chr))
  
  split(gff, gff$chr)
}

parse_temp <- function(data, r_squared_cutoff) {
  data = data %>% arrange(Dist_bp)
  linked = data %>% filter(R.2 >= r_squared_cutoff)
  unlinked = data %>% filter(R.2 < r_squared_cutoff)
  
  if (nrow(unlinked) == nrow(data)) {
    return = unlinked[1,]
  } else {
    return = linked %>% filter(R.2 >= r_squared_cutoff)
  }
  return
}

parse_chunk <- function(block, r_squared_cutoff) {
  block <- block %>% mutate(linkedSNP_count = nrow(block))
  if (block[1,]$linkedSNP_count > 1) {
    # count the number of negative/positive effects in each block
    negative <- sum(block$SNP2_effect < 0)
    positive <- sum(block$SNP2_effect > 0)
    
    # Find SNP with largest negative or positive effect
    if (positive > negative) {
      # sort in descending order
      block <- block %>% arrange(desc(SNP2_effect), Dist_bp)
      
      # get top row
      block <- block[1,]
      
    } else if (negative > positive) {
      # sort in ascending order
      block <- block %>% arrange(SNP2_effect, Dist_bp)
      
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
  } else {
    block <- block %>% mutate(linkedSNP_count = ifelse(R.2 < r_squared_cutoff, 0, 1))
  }
  block
}

tag_SNPs <- function(untagged_genes) {
  neg_genes <- untagged_genes %>% 
    filter(Effect < 0) %>% 
    arrange(Effect, P_value) %>% 
    mutate(linkedSNP_count = ifelse(linkedSNP_count <= 0, 1 ,linkedSNP_count))
  pos_genes <- untagged_genes %>% 
    filter(Effect > 0) %>% 
    arrange(desc(Effect), P_value) %>% 
    mutate(linkedSNP_count = ifelse(linkedSNP_count <= 0, 1 ,linkedSNP_count))
  negative <- sum(neg_genes$linkedSNP_count)
  positive <- sum(pos_genes$linkedSNP_count)
  
  
  if (positive > negative){
    untagged_genes <- untagged_genes %>% arrange(desc(Effect, P_value))
    tagSNP <- untagged_genes[1,]
  } else if(negative > positive){
    untagged_genes <- untagged_genes %>% arrange(Effect, P_value)
    tagSNP <- untagged_genes[1,]
  } else if(positive == negative){
    pos_max <- pos_genes[1,]
    neg_max <- neg_genes[1,]

    if(pos_max$Effect > abs(neg_max$Effect)){
      tagSNP <- pos_max
    }
    else{
      tagSNP <- neg_max
    }
  }
  tagSNP %>% mutate(linkedSNP_count=negative+positive)
}

find_genes <- function(gff, snp_df, window) {
  
  no_linked = snp_df %>% filter(linkedSNP_count == 0) %>% mutate(Position = Position1)
  one_linked = snp_df %>% filter(linkedSNP_count == 1) %>% mutate(Position = ifelse(SNP1_effect < 0 & SNP2_effect < 0, 
                                                                                    ifelse(SNP1_effect != SNP2_effect,
                                                                                           ifelse(SNP1_effect < SNP2_effect,
                                                                                                  Position1,
                                                                                                  Position2),
                                                                                           Position2),
                                                                                    ifelse(SNP1_effect > 0 & SNP2_effect > 0,
                                                                                           ifelse(SNP1_effect != SNP2_effect,
                                                                                                  ifelse(SNP1_effect > SNP2_effect,
                                                                                                         Position1,
                                                                                                         Position2),
                                                                                                  Position2),
                                                                                           ifelse(SNP1_pval != SNP2_pval,
                                                                                                  ifelse(SNP1_pval < SNP2_pval,
                                                                                                         Position1,
                                                                                                         Position2),
                                                                                                  "problem"
                                                                                           )
                                                                                    )
  )
  )
  
  more_linked = snp_df %>% filter(linkedSNP_count > 1) %>% mutate(Position = ifelse(SNP1_effect < 0 & SNP2_effect < 0,
                                                                                    ifelse(Dist_bp < window,
                                                                                           ifelse(SNP1_effect != SNP2_effect,
                                                                                                  ifelse(SNP1_effect < SNP2_effect,
                                                                                                         Position1,
                                                                                                         Position2),
                                                                                                  Position2),
                                                                                           Position2),
                                                                                    ifelse(SNP1_effect > 0 & SNP2_effect > 0,
                                                                                           ifelse(Dist_bp < window,
                                                                                                  ifelse(SNP1_effect != SNP2_effect,
                                                                                                         ifelse(SNP1_effect > SNP2_effect,
                                                                                                                Position1,
                                                                                                                Position2),
                                                                                                         Position2),
                                                                                                  Position2),
                                                                                           Position2
                                                                                    )
  )
  )
  
  problem_linked = snp_df %>% filter(linkedSNP_count == -1) %>% mutate(Position = ifelse(SNP1_effect < 0 & SNP2_effect < 0, 
                                                                                         ifelse(SNP1_effect != SNP2_effect,
                                                                                                ifelse(SNP1_effect < SNP2_effect,
                                                                                                       Position1,
                                                                                                       Position2),
                                                                                                Position2),
                                                                                         ifelse(SNP1_effect > 0 & SNP2_effect > 0,
                                                                                                ifelse(SNP1_effect != SNP2_effect,
                                                                                                       ifelse(SNP1_effect > SNP2_effect,
                                                                                                              Position1,
                                                                                                              Position2),
                                                                                                       Position2),
                                                                                                "problem")
  )
  )
  
  snp_df = rbind(no_linked, one_linked, problem_linked, more_linked) %>% 
    filter(Position != "problem") %>% mutate(Position = as.integer(Position)) %>% 
    mutate(window_start = Position-1000, window_end = Position + 1000) %>% mutate(chr = as.character(Locus1)) %>%
    mutate(Effect = ifelse(Position == Position1, SNP1_effect, SNP2_effect), P_value = ifelse(Position == Position1, SNP1_pval, SNP2_pval)) %>%
    select(chr, Position, window_start, window_end, Effect, P_value, linkedSNP_count)
  inner_join(snp_df, gff, by="chr") %>% filter((window_start >= start & window_end <= end) | (start >= window_start & end <= window_end) | 
                                                 (window_start >= start & window_start <= end) | (window_end >= start & window_end <= end)) %>% 
    mutate(Position = as.integer(Position), Effect = as.numeric(Effect), Locus1 = chr) %>%
    select(-c(chr, start, end, window_start, window_end))
}

parse_SNP <- function(all_data, LD, gff_file, window, r_squared_cutoff, num_cores) {
  
  full_gff <- read_gff(gff_file)
  
  cl <- parallel::makeCluster(num_cores, outfile="")
  registerDoParallel()
  
  all_genes = NULL
  
  # UP/DOWNSTREAM LOOP
  for (i in 1:length(LD)) {
    print(ifelse(i == 1, "upstream", "dowstream"))
    LD_stream <- LD[[i]]
    # BEGIN PROCESSING BY CHROMOSOMES LOOP
    for (name in names(LD_stream)) {
      print(name)
      if (name %in% names(full_gff)) {
        temp_data <- LD_stream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>%
          mutate(Marker2 = paste0("S", Locus1, "_", Position2)) %>% arrange(Position1)
        
        temp_data_list <- split(temp_data, temp_data$Position1)
      
        print("Parsing unlinked...")
        temp_data <- foreach(data=temp_data_list, .combine = rbind, .packages = c('dplyr', 'past')) %dopar% {
          parse_temp(data, r_squared_cutoff)
        }
        
        chr_data <- all_data %>% filter(Chr == as.integer(name))
        
        
        print("Getting effects...")
        # look up p-value and effect data for SNP1
        temp_data <- merge(temp_data, chr_data, by.x = "Marker1", by.y = "Marker") %>% 
          mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
          select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2)
        
        # look up p-value and effect data for SNP2
        temp_data <- merge(temp_data, chr_data, by.x = "Marker2", by.y = "Marker") %>%
          mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
          select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect) %>% 
          arrange(Position1)
        
        index <- c(0, cumsum(abs(diff(temp_data$Site2)) > 1))
        temp_data_list <- split(temp_data, paste(temp_data$Position1, index))
        
        print("Parsing chunks...")
        temp_data <- foreach(data=temp_data_list, .combine = rbind, .packages = c('dplyr', 'past')) %dopar% {
          parse_chunk(data, r_squared_cutoff)
        }
        
        split = 1000
        SNPS_list <- split(temp_data, rep(1:split, length.out = nrow(temp_data), each = ceiling(nrow(temp_data)/split)))
        
        # initialize genes dataframe
        SNP_genes = NULL
        
        # subset gff to only handle this chromosome
        gff <- full_gff[[as.character(name)]]
        
        # get genes in parallel
        print("Finding genes...")
        chr_genes <- foreach(SNPs=SNPS_list, .combine = rbind, .packages = c('dplyr', 'past')) %dopar%{
          find_genes(gff, SNPs, window)
        }
        
        all_genes = rbind(all_genes, chr_genes)
      }
    }
  }

  group_genes <- split(all_genes, f = all_genes$name)
  tagged_genes <- foreach(block=group_genes, .combine = rbind, .packages = c('dplyr', 'past')) %dopar% {
    tag_SNPs(block)
  }
   
  tagged_genes %>% mutate(Chromosome = Locus1, Gene = name) %>% select(Chromosome, Position, Gene, Effect, P_value)
}