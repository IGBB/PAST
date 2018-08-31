parse_pathways <- function(tagSNPs, pathways_file, gene_number_cutoff) {

  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")
  
  # UP/DOWNSTREAM LOOP
  for (i in 1:length(tagSNPs)) {
    tagSNPs_stream <- tagSNPs[[i]]

    # BEGIN PROCESSING BY CHROMOSOMES LOOP
    for (name in names(tagSNPs_stream)) {
      temp_tagSNPs <-tagSNPs_stream[[name]] 
      
      # merge pathways and temp_tagSNPs by gene
      temp_tagSNPs <- merge(pathways, temp_tagSNPs, by.x = "gene_id", by.y = "name")
      
      # filter out pathways with less than five genes associated with them
      pathways_filtered <- temp_tagSNPs %>% group_by(pathway_id) %>% dplyr::summarise(count = n()) %>% filter(count >= gene_number_cutoff)
      temp_tagSNPs <- temp_tagSNPs %>% filter(pathway_id %in% pathways_filtered$pathway_id)
      
      # store pathway information for chromosomes
      tagSNPs_stream[[name]] <- temp_tagSNPs
      
    } # END CHROMOSOMES LOOP
    
    # store pathway information for up/dowstream
    tagSNPs[[i]] <- tagSNPs_stream
    
  } # END UP/DOWNSTREAM LOOP
  
# return updated tagSNPs with pathway information
tagSNPs
}