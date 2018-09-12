parse_pathways <- function(genes, pathways_file, gene_number_cutoff) {

  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")

  # merge pathways and temp_tagSNPs by gene
  genes <- merge(pathways, genes, by.x = "gene_id", by.y = "name")
  
  # filter out pathways with less than five genes associated with them
  pathways_filtered <- genes %>% group_by(pathway_id) %>% dplyr::summarise(count = n()) %>% filter(count >= gene_number_cutoff)
  
  # return filtered pathways
  genes %>% filter(pathway_id %in% pathways_filtered$pathway_id)

}