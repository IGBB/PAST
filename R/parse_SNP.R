library(dplyr)
read_gff <- function(gff_file) {
  read.table(gff_file, comment.char = "#", col.names = 
                      c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
    filter(type == "gene") %>% 
    select(chr, attributes, start, end) %>% 
    mutate(name = strsplit(as.character(strsplit(as.character(attributes), ";Name=")[[1]][2]), ";")[[1]][1], attributes=NULL)
}

find_gene <- function(gff, snp_df) {
  window_start <- snp_df$Position1-1000
  window_end <- snp_df$Position1+1000
  snp_chr <- snp_df$Locus1
  snp_marker <- snp_df$Marker1
  
  gene <- gff %>% filter(chr == snp_chr) %>% 
    filter((window_start >= start & window_end <= end) | (start >= window_start & end <= window_end) | 
             (window_start >= start & window_start <= end) | (window_end >= start & window_end <= end))
  
  if (nrow(gene) > 0) { 
    gene <- gene %>% mutate(Marker1 = snp_marker)
    merge(snp_df, gene, by = "Marker1") #%>% mutate(chr = NULL, start = NULL, end = NULL)
  } else {
    mutate(snp_df, name="NA")
  }
}

parse_SNP <- function(all_data, LD, gff_file) {
  # This line will eventually be the line below, but for now, we can hardcode the GFF file.
  # gff <- read_gff(gff_file)
  gff <- read_gff("example/zea_mays.protein_coding.gff.gz")

  # get upstream and downstream
  LD_upstream <- LD[[1]]
  LD_downstream <- LD[[2]]

  # BEGIN UPSTREAM LOOP
  for (name in names(LD_upstream)) {
    temp_data <- LD_upstream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>%
      mutate(Marker2 = paste0("S", Locus1, "_", Position2))

    # filter based on SNPs in stats/effects
    temp_data <- temp_data %>% filter(temp_data$Position1 %in% all_data$Pos)

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
          blocks[[block_name]] <- find_gene(gff, block)
        } else if (negative > positive) {
          # sort in ascending order
          block <- block %>% arrange(SNP2_effect)

          # get top row
          block <- block[1,]

          # store gene
          blocks[[block_name]] <- find_gene(gff, block)
        } else if (negative == positive) {
          if(block$SNP1_effect[[1]] > 0){
            # sort in descending order
            block <- block %>% arrange(desc(SNP2_effect))
            
            # get top row
            block <- block[1,]
            
            # store gene
            blocks[[block_name]] <- find_gene(gff, block)
          }
          else{
            # sort in ascending order
            block <- block %>% arrange(SNP2_effect)
            
            # get top row
            block <- block[1,]
            
            # store gene
            blocks[[block_name]] <- find_gene(gff, block)
          }
        }
      } # END BLOCKS LOOP
      
      # rejoin blocks and overwrite block_SNPs
      # you'll need to run the whole for block in blocks loop
      # to get blocks set up correctly
      blocks <- rbind.fill(blocks)
      block_genes <- rbind(block_genes, blocks)
      
    } # END POSITIONS SNPS LOOP

    # get all the single SNPs in a data frame by themselves
    single_SNPs <- chr_linked %>% group_by(Position1) %>% dplyr::summarise(count = n()) %>% filter(count == 1)
    single_SNPs <- chr_linked %>% filter(Position1 %in% single_SNPs$Position1)
    
    # we can condense our code a little by using the variable in the loop directly
    # instead of assigning it
    # positions <- single_SNPs$Position1
    
    for(pos in single_SNPs$Position1){
      row <- single_SNPs %>% filter(Position1 == pos)
      single_SNPs[[pos]] <- find_gene(gff, row)
    }

    # get genes here

    # catch everything that didn't make it through the filter
    chr_unlinked <- temp_data %>% filter(R.2 >= 0.8) %>% arrange(Position1, Dist_bp)

    # find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]

    # get genes here

    ## Last two steps are below this line. All your stuff should be above this.

    # combine linked and unlinked and arrange Position1 in ascending order
    temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)

    # store modified and filtered data
    LD_upstream[[name]] <- temp_data

  } # END UPSTREAM LOOP  

  # BEGIN DOWNSTREAM LOOP
  for (name in names(LD_downstream)) {
    temp_data <- LD_upstream[[name]] %>% mutate(Marker1 = paste0("S", Locus1, "_", Position1)) %>%
      mutate(Marker2 = paste0("S", Locus1, "_", Position2))

    # filter based on SNPs in stats/effects
    temp_data <- temp_data %>% filter(temp_data$Position1 %in% all_data$Pos)

    ## BEGIN PROCESSING BLOCKS OF LINKED SNPS

    # retrieve linked SNPs
    chr_linked <- temp_data %>% arrange(Position1) %>% filter(R.2 >= 0.8)

    # get SNPs that need to be parsed by block
    block_SNPs <- chr_linked %>% group_by(Marker1) %>% dplyr::summarise(count = n()) %>% filter(count > 1)
    block_SNPs <- chr_linked %>% filter(Marker1 %in% block_SNPs$Marker1)

    # find all unique positions in block_SNPs
    positions <- block_SNPs %>% group_by(Position1) %>% dplyr::summarise(count = n())

    # begin block SNP loops
    for(pos in positions$Position1) {
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
        block <- merge(block, all_data, by.x = "Marker2", by.y = "Marker") %>% 
        mutate(SNP1_pval = p, SNP1_effect = Effect.x) %>% 
        select(Locus1, Position1, Position2, Site1, Site2, Dist_bp, R.2, Marker1, SNP1_pval, SNP1_effect, Marker2, linkedSNP_count)

        # look up p-value and effect data for SNP2
        block <- merge(block, all_data, by.x = "Marker1", by.y = "Marker") %>%
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
          blocks[[block_name]] <- find_gene(gff, block)
        } else if (negative > positive) {
          # sort in ascending order
          block <- block %>% arrange(SNP2_effect)

          # get top row
          block <- block[1,]

          # store gene
          blocks[[block_name]] <- find_gene(gff, block)
        } else if (negative == positive) {
          if (block$SNP1_effect >0){
            # sort in descending order
            block <- block %>% arrange(desc(SNP2_effect))
            
            # get top row
            block <- block[1,]
            
            # store gene
            blocks[[block_name]] <- find_gene(gff, block)
          }
          else{
            # sort in ascending order
            block <- block %>% arrange(SNP2_effect)
            
            # get top row
            block <- block[1,]
            
            # store gene
            blocks[[block_name]] <- find_gene(gff, block)
          }
        }
      } # END BLOCKS LOOP
    } # END POSITIONS SNPS LOOP

    # get all the single SNPs in a data frame by themselves
    single_SNPs <- chr_linked %>% group_by(Position1) %>% dplyr::summarise(count = n()) %>% filter(count == 1)

    # get genes here

    # catch everything that didn't make it through the filter
    chr_unlinked <- temp_data %>% filter(R.2 >= 0.8) %>% arrange(Position1, Dist_bp)

    # find the first instance of every unlinked SNP
    chr_unlinked <- chr_unlinked[match(unique(chr_unlinked$Position1), chr_unlinked$Position1),]

    # get genes here

    ## Last two steps are below this line. All your stuff should be above this.

    # combine linked and unlinked and arrange Position1 in ascending order
    temp_data <- rbind(chr_linked, chr_unlinked) %>% arrange(Position1)

    # store modified and filtered data
    LD_downstream[[name]] <- temp_data

  } # END DOWNSTREAM LOOP

  # return modified lists
  list(LD_upstream, LD_downstream)
}