#Get Factors
get_factors <- function(gene_ranks) {
  factors = vector("logical", length(gene_ranks))
  factors[1] = gene_ranks[1] - 1
  for (i in 2:length(gene_ranks)) {
    factors[i] = gene_ranks[i] - gene_ranks[i-1] - 1 + factors[i-1]
  }
  factors
}

find_significant_pathways <- function(genes, pathways_file, gene_number_cutoff, mode) {
  
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
    if (mode == "resistant") {
      temp_data <- temp_data %>% arrange(Effect) %>% mutate(rank = row_number())
    } else if (mode == "susceptible") {
      temp_data <- temp_data %>% arrange(desc(Effect)) %>% mutate(rank = row_number())
    } else {
      stop("Incorrect mode.")
    }
    
    column_observations = data.frame()
    
    foreach(pathway = iter(pathways_unique$pathway_id, by='row')) %do% {
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
        # this function needs to return a dataframe with
        #Phits_Pmisses <- get_Phits_Pmisses(Phits, Pmisses)
        
        # I don't think you need a function for max
        # you can just sort descending (I think) and take the
        # top row
        #max <- most negative Phits_Pmisses
        
        column_observations <- rbind(column_observations, sort(factors, decreasing = TRUE)[1])
      } else {
        column_observations <- rbind(column_observations, NA)
      }
    }
    pathways_unique <- cbind(pathways_unique, column_observations)
  }
  colnames(pathways_unique) <- c("Pathway", "ES_Observed", 1:1000)
  colnames(pathways_unique)[3:1002] <- paste0("ES", colnames(pathways_unique)[3:1002])
  pathways_unique <- pathways_unique %>% filter(!is.na(ES_Observed))
  pathways_unique <- pathways_unique %>% mutate(permutation_mean = apply(pathways_unique[,3:1002], 1, mean))
  pathways_unique <- pathways_unique %>% mutate(permutation_standard_deviation = apply(pathways_unique[,3:1002], 1, sd))
  pathways_unique <- pathways_unique %>% mutate(NES_Observed = (ES_Observed-permutation_mean)/permutation_standard_deviation)
  for (i in 1:1000) {
    pathways_unique[i+2] = (pathways_unique[i+2] - pathways_unique[7]) / pathways_unique[8]
  }
  pathways_unique <- pathways_unique %>% mutate(pvalue = 1-pnorm(NES_Observed))
  
  FDR_denominators<-rep(NA,length(pathways_unique$NES_Observed))
  for(i in 1:length(pathways_unique$NES_Observed)){
    FDR_denominators[i]<-(sum(pathways_unique$NES_Observed>=pathways_unique$NES_Observed[i]))/length(pathways_unique$NES_Observed)
  }
  
  FDR_numerators <- rep(NA,length(pathways_unique$NES_Observed))
  effects =
  for(i in 1:length(pathways_unique$NES_Observed)) {
    FDR_numerators[i] <- (sum(pathways_unique$NES_Observed>=pathways_unique$ES_Observed[i]))/(ncol(pathways_unique[3:1002])*length(pathways_unique$NES_Observed))
  }
  
  pathways_unique <- mutate(pathways_unique, FDR = FDR_numerators/FDR_denominators) %>% arrange(pvalue)
  pathways_unique <- pathways_unique %>% select(Pathway, ES_Observed, NES_Observed, pvalue, FDR) %>%
    mutate(qvalue = qvalue(pathways_unique$pvalue, lambda=0, fdr.level = 0.05)$qvalues)
  
  pathways_significant <- pathways_unique %>% filter(pvalue <= 0.2) %>% select(Pathway, NES_Observed) %>% mutate(NESrank = row_number())
  temp_data <- as.data.frame(cbind(effects[,1], effects[,2]))
  colnames(temp_data) <- c("Gene", "Effect")
  if (mode == "resistant") {
    temp_data <- temp_data %>% arrange(Effect) %>% mutate(rank = row_number())
  } else if (mode == "susceptible") {
    temp_data <- temp_data %>% arrange(desc(Effect)) %>% mutate(rank = row_number())
  } else {
    stop("Incorrect mode.")
  }
  
  rugplots_data = NULL
  foreach(pathway = iter(pathways_significant$pathway_id, by='row')) %do% {
    genes_in_pathway <- filter(pathways, pathway_id == pathway)
    
    ## get ranks and effects and sort by rank
    genes_in_pathway <- merge(genes_in_pathway, temp_data, by.x = "gene_id", by.y = "Gene") %>% arrange(rank)
    
    # get factors using rank
    genes_in_pathway$factors <- get_factors(genes_in_pathway$rank) 
    
    # you need a get_Pmisses function
    #genes_in_pathway$Pmisses <- get_Pmisses(temp_data, genes_in_pathway, factors)
    
    # since the NR is just a summation, do that at the
    # beginning of this function. R can sum a column easily.
    # you need to create this function
    #genes_in_pathway$Phits <- get_Phits(genes_in_pathway$effect)
    
    # create this function too
    # this function needs to return a dataframe with
    #genes_in_pathway$phit_pmiss <- get_Phits_Pmisses(Phits, Pmisses)
    
    # append rows with NESrank
    genes_in_pathway<- merge(genes_in_pathway, pathways_significant, by.x="pathway_id", by.y="Pathway")
    rugplots_data <- rbind(rugplots_data, genes_in_pathway)
  }
  
  rugplots_data <- rugplots_data %>% select(pathway_id, NESrank, gene_id, rank, phit_pmiss)
}
