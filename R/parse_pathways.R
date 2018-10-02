find_significant_pathways <- function(genes, pathways_file, gene_number_cutoff) {

  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")

  # merge pathways and temp_tagSNPs by gene
  genes <- merge(pathways, genes, by.x = "gene_id", by.y = "Gene") %>% mutate(Gene = gene_id, gene_id = NULL)
  
  # filter out pathways with less than five genes associated with them
  pathways_filtered <- genes %>% group_by(pathway_id) %>% dplyr::summarise(count = n()) %>% filter(count >= gene_number_cutoff)
  
  # return filtered pathways
  genes_with_pathways <- genes %>% filter(pathway_id %in% pathways_filtered$pathway_id)
  
  # sample to create 1000 random distributions
  effects <- genes %>% select(Gene, Effect)
  effects <- cbind(effects, sapply(1:1000, function(i) sample(effects$Effect)))
  
  pathways_unique <- unique(select(pathways, pathway_id))
  pathways_unique[] <- lapply(pathways_unique, as.character)
  
  # in the original code, Juliet gets the gene IDs in their own dataset
  # we can always access them with effects[,1] so we don't have to do that
  for (i in 2:1003) {
    temp_data <- as.data.frame(cbind(effects[,1], effects[,i]))
    colnames(temp_data) <- c("Gene", "Effect")
    temp_data <- temp_data %>% arrange(Effect) %>% mutate(Rank = row_number())
    
    foreach(pathway = iter(pathways_unique, by='row')) %do% {
      print(pathway)
      genes_in_pathway <- filter(pathways, pathway_id == pathway)
      
      ## getSranks(genes_in_pathway, temp_data)
      # this particular part of the code gets
      # the rank of the genes in the pathway
      # from the rank of the genes in the whole
      # dataset. You can just merge the two
      # datasets. This probably doesn't need a
      # whole function in R.
      genes_in_pathway <- <YOUR CODE HERE>
      
      ## getSeffects(genes_in_pathway, temp_data)
      # this part of the code gets the effects of
      # the genes in the dataset from the whole
      # dataset. If you did the merge above, you
      # can just keep the effects data when you
      # merge instead of doing anything extra
      genes_in_pathway <- <YOUR CODE HERE>
      
      ## check cutoff
      if <YOUR CODE HERE> {
        
        # you need to write a get_factors function that takes
        # the rank column as input
        factors <- get_factors(genes_in_pathway$rank)
        
        # you need a get_Pmisses function
        Pmisses <- get_Pmisses(temp_data, genes_in_pathway, factors)
        
        # since the NR is just a summation, do that at the
        # beginning of this funciton. R can sum a column easily.
        # you need to create this function
        Phits <- get_Phits(genes_in_pathway$effect)
        
        # create this function too
        Phits_Pmisses <- get_Phits_Pmisses(Phits, Pmisses)
        
        # I don't think you need a function for max
        # you can just sort descending (I think) and take the
        # top row
        max <- #most negative Phits_Pmisses
          
        # I'll work on the output later
        
        
      }
    }
    
  }
  
  
  
  # 11.R code goes here
}