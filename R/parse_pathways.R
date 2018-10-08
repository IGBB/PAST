#Get Factors
get_factors <- function(gene_ranks) {
  factors = vector("logical", length(gene_ranks))
  factors[1] = gene_ranks[1] - 1
  for (i in 2:length(gene_ranks)) {
    factors[i] = gene_ranks[i] - gene_ranks[i-1] - 1 + factors[i-1]
  }
  factors
}

find_significant_pathways <- function(genes, pathways_file, gene_number_cutoff) {

  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")
  
  # sample to create 1000 random distributions
  effects <- genes %>% select(Gene, Effect)
  effects <- cbind(effects, sapply(1:1000, function(i) sample(effects$Effect)))
  effects$Gene <- as.character(effects$Gene)
  
  pathways_unique <- unique(select(pathways, pathway_id))
  pathways_unique[] <- lapply(pathways_unique, as.character)
  
  # process sample columns
  for (i in 2:1002) {
    print(i)
    temp_data <- as.data.frame(cbind(effects[,1], effects[,i]))
    colnames(temp_data) <- c("Gene", "Effect")
    temp_data <- temp_data %>% arrange(Effect) %>% mutate(rank = row_number())
    
    column_observations = data.frame()
    
    foreach(pathway = iter(pathways_unique$pathway_id, by='row')) %do% {
      print(pathway)
      genes_in_pathway <- filter(pathways, pathway_id == pathway)
      
      ## get ranks and effects and sort by rank
      genes_in_pathway <- merge(genes_in_pathway, temp_data, by.x = "gene_id", by.y = "Gene") %>% arrange(rank)
      
      ## check cutoff
      if (nrow(genes_in_pathway) > gene_number_cutoff) {
        
        # get factors using rank
        factors <- get_factors(genes_in_pathway$rank) 
        
        # you need a get_Pmisses function
        #Pmisses <- get_Pmisses(temp_data, genes_in_pathway, factors)
        
        # since the NR is just a summation, do that at the
        # beginning of this function. R can sum a column easily.
        # you need to create this function
        #Phits <- get_Phits(genes_in_pathway$effect)
        
        # create this function too
        #Phits_Pmisses <- get_Phits_Pmisses(Phits, Pmisses)
        
        # I don't think you need a function for max
        # you can just sort descending (I think) and take the
        # top row
        #max <- most negative Phits_Pmisses
        
        column_observations <- rbind(column_observations, max)
      } else {
        column_observations <- rbind(column_observations, NA)
      }
    }
    pathways_unique <- cbind(pathways_unique, column_observations)
  }
  colnames(pathways_unique) <- c("Pathway", "ES_Observed", 1:1000)
  
  # 11.R code goes here
}
