read_gff <- function(gff_file) {

  read.table(gff_file, comment.char = "#", col.names = 
                      c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))  %>% 
    select(chr, attributes, start, end) %>% 
    mutate(tempcol = stringr::str_split(attributes, ';Name=')) %>%
    rowwise() %>%
    mutate(name = gsub("ID=", "", unlist(tempcol)[1]), tempcol = NULL, attributes = NULL)
}

find_gene <- function(gff, snp_df, window) {
  problem = FALSE
  snp_df$Dist_bp <- as.numeric(levels(snp_df$Dist_bp))[snp_df$Dist_bp]
  if (snp_df$linkedSNP_count == 0) {
    position = snp_df$Position1
  } else if (snp_df$linkedSNP_count == -1 | snp_df$linkedSNP_count == 1){
      if (snp_df$SNP1_effect != snp_df$SNP2_effect) {
        if (snp_df$SNP1_effect < 0 & snp_df$SNP2_effect < 0) {
          if (snp_df$SNP1_effect < snp_df$SNP2_effect) {
            position = snp_df$Position1
          } else {
            position = snp_df$Position2
          }
        } else if (snp_df$SNP1_effect > 0 & snp_df$SNP2_effect > 0) {
          if (snp_df$SNP1_effect > snp_df$SNP2_effect) {
            position = snp_df$Position1
          } else {
            position = snp_df$Position2
          }
        } else {
          if (snp_df$linkedSNP_count == -1) {
            print("ERROR: effects have opposite signs")
          } else {
              if (snp_df$SNP1_pval != snp_df$SNP2_pval) {
                if (snp_df$SNP1_pval < snp_df$SNP2_pval) {
                  position = snp_df$Position1
                } else {
                  position = snp_df$Position2
                }
              } else {
                print("Single-linked SNP error: same p-value")
                problem = TRUE
                position = snp_df$Position2
              }
          }
        }
      } else {
        position = snp_df$Position2
      }
  } else {
    if (snp_df$Dist_bp <= window) {
      if (snp_df$SNP1_effect == snp_df$SNP2_effect) {
        position = snp_df$Position2
      } else if (snp_df$SNP1_effect < 0 & snp_df$SNP2_effect < 0) {
        if (snp_df$SNP1_effect < snp_df$SNP2_effect) {
          position = snp_df$Position1
        } else {
          position = snp_df$Position2
        }
      } else if (snp_df$SNP1_effect > 0 & snp_df$SNP2_effect > 0) {
        if (snp_df$SNP1_effect > snp_df$SNP2_effect) {
          position = snp_df$Position1
        } else {
          position = snp_df$Position2
        }
      } else {
        if (snp_df$SNP1_pval < snp_df$SNP2_pval) {
          position = snp_df$Position1
        } else {
          position = snp_df$Position2
        }
      }
    } else {
      position = snp_df$Position2
    }
  }
  window_start <- position-window
  window_end <- position+window
  snp_chr <- snp_df$Locus1
  snp_marker <- snp_df$Marker1
  
  gene <- gff %>% filter(chr == snp_chr) %>% 
    filter((window_start >= start & window_end <= end) | (start >= window_start & end <= window_end) | 
             (window_start >= start & window_start <= end) | (window_end >= start & window_end <= end))
  
  if (nrow(gene) > 0 & problem == FALSE) { 
    gene <- gene %>% mutate(Marker1 = snp_marker)
    merge(snp_df, gene, by = "Marker1") %>% mutate(chr = NULL, start = NULL, end = NULL)
  } else {
    mutate(snp_df, name="NA") %>% 
      select(Marker1, Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count, name)
  }
}

parse_SNP <- function(all_data, LD, gff_file, window) {
  # This line will eventually be the line below, but for now, we can hardcode the GFF file.
  # gff <- read_gff(gff_file)
  gff <- read_gff("example/example.gff.gz")
  
  # UP/DOWNSTREAM LOOP
  for (i in 1:length(LD)) {
    LD_stream <- LD[[i]]
    names(LD_stream)
  
    # BEGIN PROCESSING BY CHROMOSOMES LOOP
    for (name in names(LD_stream)) {
      temp_data <- LD_stream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>%
        mutate(Marker2 = paste0("S", Locus1, "_", Position2))
      
      # temporary line to subset data for testing
      temp_data <- head(temp_data, 5000)
  
      # filter based on SNPs in stats/effects
      # temp_data <- temp_data %>% filter(temp_data$Position1 %in% all_data$Pos)
  
      # retrieve linked SNPs
      chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= 0.8)
  
      # get SNPs that need to be parsed by block
      block_SNPs <- chr_linked %>% group_by(Marker1) %>% dplyr::summarise(count = n()) %>% filter(count > 1)
      block_SNPs <- chr_linked %>% filter(Marker1 %in% block_SNPs$Marker1)
  
      # find all unique positions in block_SNPs
      positions <- block_SNPs %>% group_by(Position1) %>% dplyr::summarise(count = n())
  
      # initialize block genes
      block_genes = NULL
      
      # begin block SNP loops
      for(pos in positions$Position1) {
        print(pos)
        # select all SNPs at pos
        positions_block_SNPs <- block_SNPs %>% filter(Position1 == pos)
  
        # find indices where the difference in Site2 > 1 and split into blocks
        index <- c(0, cumsum(abs(diff(positions_block_SNPs$Site2)) > 1))
        blocks <- split(positions_block_SNPs, index)
  
        # process blocks
        for (block_name in names(blocks)) {
          
          # get current block
          block <- blocks[[block_name]]
          # add linkedSNP_count
          block <- block %>% mutate(linkedSNP_count = nrow(block))
  
          # look up p-value and effect data for SNP1
          block <- merge(block, all_data, by.x = "Marker1", by.y = "Marker") %>% 
            mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
            select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)
  
          # look up p-value and effect data for SNP2
          block <- merge(block, all_data, by.x = "Marker2", by.y = "Marker") %>%
            mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
            select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count)
  
          # count the number of negative/positive effects in each block
          negative <- sum(block$SNP2_effect < 0)
          positive <- sum(block$SNP2_effect > 0)
  
          # Find SNP with largest negative or positive effect
          if (positive > negative) {                 
          # sort in descending order
            block <- block %>% arrange(desc(SNP2_effect))
  
            # get top row
            block <- block[1,]
            # store gene
            blocks[[block_name]] <- find_gene(gff, block, window)
          } else if (negative > positive) {
            # sort in ascending order
            block <- block %>% arrange(SNP2_effect)
  
            # get top row
            block <- block[1,]
  
            # store gene
            blocks[[block_name]] <- find_gene(gff, block, window)
          } else if (negative == positive) {
            if(block$SNP1_effect[[1]] > 0){
              # sort in descending order
              block <- block %>% arrange(desc(SNP2_effect))
              
              # get top row
              block <- block[1,] %>% mutate(linkedSNP_count = -1)
              
              # store gene
              blocks[[block_name]] <- find_gene(gff, block, window)
            }
            else{
              # sort in ascending order
              block <- block %>% arrange(SNP2_effect)
              
              # get top row
              block <- block[1,] %>% mutate(linkedSNP_count = -1)
              
              # store gene
              blocks[[block_name]] <- find_gene(gff, block, window)
            }
          }
        } # END BLOCKS LOOP
        
        # rejoin blocks and overwrite block_SNPs
        # you'll need to run the whole for block in blocks loop
        # to get blocks set up correctly
        blocks <- rbind.fill(blocks)
        block_genes <- rbind(block_genes, blocks)
        
      } # END POSITIONS SNPS LOOP
      
      block_genes <- block_genes %>% filter(name != "NA")
      
      # get all the single SNPs in a data frame by themselves and get SNP data
      single_SNP_genes = NULL
      single_SNPs <- chr_linked %>% group_by(Marker1) %>% dplyr::summarise(count = n()) %>% filter(count == 1)
      single_SNPs <- chr_linked %>% filter(Marker1 %in% single_SNPs$Marker1) %>% mutate(linkedSNP_count = 1)
      single_SNPs <- merge(single_SNPs, all_data, by.x = "Marker1", by.y = "Marker") %>% 
        mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)
      single_SNPs <- merge(single_SNPs, all_data, by.x = "Marker2", by.y = "Marker") %>%
        mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count)
      
      # get genes for single linked SNPs and filter NA genes
      for(pos in single_SNPs$Position1){
        print(pos)
        row <- single_SNPs %>% filter(Position1 == pos)
        gene <- find_gene(gff, row, window)
        single_SNP_genes <- rbind(single_SNP_genes, gene)
      }
      
      single_SNP_genes <- single_SNP_genes %>% filter(name != "NA")
      
      test <- single_SNP_genes %>% group_by(name) %>% dplyr::summarise(count = n())
  
      # catch everything that didn't make it through the filter
      chr_unlinked <- temp_data %>% filter(R.2 < 0.8) %>% arrange(Position1, Dist_bp)
  
      # find the first instance of every unlinked SNP
      unlinked_SNPs <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]

      # get genes for unlinked SNPs and filter NA genes
      unlinked_SNP_genes = NULL
      unlinked_SNPs <- unlinked_SNPs %>% mutate(linkedSNP_count = 0)
      unlinked_SNPs <- merge(unlinked_SNPs, all_data, by.x = "Marker1", by.y = "Marker") %>% 
        mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)
      unlinked_SNPs <- merge(unlinked_SNPs, all_data, by.x = "Marker2", by.y = "Marker") %>%
        mutate(SNP2_pval = p, SNP2_effect = Effect.x) %>%
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, SNP2_pval, SNP2_effect, linkedSNP_count)
      
      for(pos in unlinked_SNPs$Position1){
        print(pos)
        row <- unlinked_SNPs %>% filter(Position1 == pos)
        gene <- find_gene(gff, row, window)
        unlinked_SNP_genes <- rbind(unlinked_SNP_genes, gene)
      }
      
      unlinked_SNP_genes <- unlinked_SNP_genes %>% filter(name != "NA")
  
      ## Last two steps are below this line. All your stuff should be above this.
  
      # combine linked and unlinked and arrange Position1 in ascending order
      temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)
  
      # store modified and filtered data
      LD_stream[[name]] <- temp_data
  
    } # END UP/DOWNSTREAM LOOP
    
    # set up/downstream to modified data
    LD[[i]] <- LD_stream
    
  }
  # Gene sorting loop
  all_genes <- rbind(block_genes, single_SNP_genes, unlinked_SNP_genes)
  all_genes <- arrange(all_genes, desc(name))
  all_genes <- all_genes %>% mutate(linkedSNP_count <- ifelse(linkedSNP_count <= 0, 1 ,linkedSNP_count))
  group_genes <- split(all_genes, f = all_genes$name)
  for ( name_block in names(group_genes)){
    single <- group_genes[[name_block]]
    neg_genes <- single %>% filter(SNP2_effect < 0)
    pos_genes <- single %>% filter(SNP2_effect > 0)
    negative <- sum(neg_genes$linkedSNP_count)
    positive <- sum(pos_genes$linkedSNP_count)
    
    if (positive > negative){
      single <- single %>% arrange(desc(SNP2_effect))
      tagSNP <- single[1,]
  } else if(negative > positive){
    single <- single %>% arrange(SNP2_effect)
    tagSNP <- single[1,]
  } else if(positive == negative){
    if (single$SNP1_effect[[1]] > 0){
      single <- single %>% arrange(desc(SNP2_effect))
      tagSNP <- single[1,]
    } 
    else{
      single <- single %>% arrange(SNP2_effect)
      tagSNP <- single[1,]
    }
      
  }
  }
  
  # return LD
  LD
}