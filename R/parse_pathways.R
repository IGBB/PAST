get_factors <- function(gene_ranks) {
  factors = vector("logical", length(gene_ranks))
  factors[1] = gene_ranks[1] - 1
  for (i in 2:length(gene_ranks)) {
    factors[i] = gene_ranks[i] - gene_ranks[i-1] - 1 + factors[i-1]
  }
  factors
}

get_pmisses <- function(temp_data, genes_in_pathway, factors) {
  pmisses = factors/(nrow(temp_data) - nrow(genes_in_pathway))
  pmisses
}

get_phits <- function(gene_effects) {
  NR = sum(abs(gene_effects))
  
  absolute = abs(gene_effects)
  phits = vector("logical", length(gene_effects))
  
  phits[1] = (absolute[1] / NR)
  for (i in 2:length(gene_effects)) {
    phits[i] = absolute[i]/NR + phits[i-1]
  }
  phits
}

get_phits_pmisses <- function(phits, pmisses) {
  phits_pmisses = phits - pmisses
  phits_pmisses
}

find_max <- function(phits_pmisses){
 phits_pmisses = sort(phits_pmisses, decreasing = TRUE)
 max = phits_pmisses[[1]]
 max
}

find_significant_pathways <- function(genes, pathways_file, gene_number_cutoff, mode, sample_size, num_cores) {
  
  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")
  
  # sample to create 1000 random distributions
  effects <- genes %>% select(Gene, Effect)
  effects <- cbind(effects, sapply(1:sample_size, function(i) sample(effects$Effect)))
  
  #effects <- read.table("/home/thrash/Dropbox/maize/effperms/Aflaenv1_effperms.txt", sep="\t", header=TRUE)
  #effects <- read.table("/home/thrash/Work/marilyn/past-test/aflatoxin/25404.effperms.txt", sep="\t", header=TRUE)
  results <- read.table("/home/thrash/Dropbox/maize/effperms/Aflaenv1v4_nop_ESq.txt", sep="\t", header=TRUE) %>%
    mutate(Pathway = as.character(PWid), PWid=NULL, NES_Observed = NESobs, NESobs = NULL, pvalue = pval, pval = NULL, qvalue = qobj.qvalues, qobj.qvalues=NULL) %>%
    select(Pathway, NES_Observed, FDR, pvalue, qvalue)

  pathways_unique <- unique(select(pathways, pathway_id))
  pathways_unique[] <- lapply(pathways_unique, as.character)
  
  # process sample columns
  cl <- parallel::makeCluster(num_cores, outfile="")
  registerDoParallel()
  column_observations <- foreach(i=iter(2:(sample_size+2)), .combine = cbind, .packages = c('dplyr', 'past', 'foreach','iterators')) %dopar%{
    temp_data <- data.frame(matrix("NA", ncol = 2, nrow = nrow(effects)))
    colnames(temp_data) <- c("Gene", "Effect")
    temp_data$Effect = effects[,i]
    temp_data$Gene = effects[,1]
    
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
      genes_in_pathway <- merge(genes_in_pathway, temp_data, by.x = "gene_id", by.y = "Gene") %>% arrange(rank) %>% unique()
      
      ## check cutoff
      if (nrow(genes_in_pathway) >= gene_number_cutoff) {
        
        # get factors using rank
        print("factors")
        factors <- get_factors(genes_in_pathway$rank) 
        
        # get pmisses
        print("pmisses")
        pmisses <- get_pmisses(temp_data, genes_in_pathway, factors)
        
        # get phits
        print("phits")
        phits <- get_phits(genes_in_pathway$Effect)
        
        # get phits-pmisses
        print("phits-pmisses")
        phits_pmisses <- get_phits_pmisses(phits, pmisses)
        find_max(phits_pmisses)
        # store max phit_pmisses
        column_observations <- rbind(column_observations, find_max(phits_pmisses))
      } else {
        column_observations <- rbind(column_observations, NA)
      }
    }
    column_observations
  }
  
  stopCluster(cl)

  pathways_unique <- cbind(pathways_unique, column_observations)
  colnames(pathways_unique) <- c("Pathway", "ES_Observed", 1:sample_size)
  colnames(pathways_unique)[3:(sample_size+2)] <- paste0("ES", colnames(pathways_unique)[3:(sample_size+2)])
  pathways_unique <- pathways_unique %>% filter(!is.na(ES_Observed))

  pathways_unique <- pathways_unique %>% mutate(permutation_mean = apply(pathways_unique[,3:(sample_size+2)], 1, mean))
  pathways_unique <- pathways_unique %>% mutate(permutation_standard_deviation = apply(pathways_unique[,3:(sample_size+2)], 1, sd))
  pathways_unique <- pathways_unique %>% mutate(NES_Observed = (ES_Observed-permutation_mean)/permutation_standard_deviation)
  for (i in 1:sample_size) {
    pathways_unique[i+2] = (pathways_unique[i+2] - pathways_unique[sample_size+3]) / pathways_unique[sample_size+4]
  }
  pathways_unique <- pathways_unique %>% mutate(pvalue = 1-pnorm(NES_Observed))
  
  FDR_denominators<-rep(NA,length(pathways_unique$NES_Observed))
  for(i in 1:length(pathways_unique$NES_Observed)){
    FDR_denominators[i]<-(sum(pathways_unique$NES_Observed>=pathways_unique$NES_Observed[i]))/length(pathways_unique$NES_Observed)
  }
  
  FDR_numerators <- rep(NA,length(pathways_unique$NES_Observed))
  for(i in 1:length(pathways_unique$NES_Observed)) {
    FDR_numerators[i] <- (sum(pathways_unique[,3:(sample_size+2)]>=pathways_unique$NES_Observed[i]))/(ncol(pathways_unique[3:(sample_size+2)])*length(pathways_unique$NES_Observed))
  }

  pathways_unique <- mutate(pathways_unique, FDR = FDR_numerators/FDR_denominators) %>% arrange(pvalue)
  pathways_unique <- pathways_unique %>% select(Pathway, ES_Observed, NES_Observed, pvalue, FDR) %>%
    mutate(qvalue = qvalue(pathways_unique$pvalue, lambda=0, fdr.level = 0.05)$qvalues)

  pathways_significant <- pathways_unique %>% select(Pathway, NES_Observed, FDR, pvalue, qvalue) %>% mutate(NESrank = row_number())
  
  temp_data <- data.frame(matrix("NA", ncol = 2, nrow = nrow(effects)))
  colnames(temp_data) <- c("Gene", "Effect")
  temp_data$Effect = effects[,2]
  temp_data$Gene = effects[,1]
  colnames(temp_data) <- c("Gene", "Effect")
  if (mode == "resistant") {
    temp_data <- temp_data %>% arrange(Effect) %>% mutate(rank = row_number())
  } else if (mode == "susceptible") {
    temp_data <- temp_data %>% arrange(desc(Effect)) %>% mutate(rank = row_number())
  } else {
    stop("Incorrect mode.")
  }

  rugplots_data = NULL
  rugplots_data <- foreach(pathway = iter(pathways_significant$Pathway, by='row'), .combine = rbind) %do% {
    genes_in_pathway <- filter(pathways, pathway_id == pathway)

    ## get ranks and effects and sort by rank
    genes_in_pathway <- merge(genes_in_pathway, temp_data, by.x = "gene_id", by.y = "Gene") %>% arrange(rank)

    # get factors using rank
    genes_in_pathway$factors <- get_factors(genes_in_pathway$rank)

    # get pmisses
    genes_in_pathway$pmisses <- get_pmisses(temp_data, genes_in_pathway, genes_in_pathway$factors)

    # get phits
    genes_in_pathway$phits <- get_phits(genes_in_pathway$Effect)

    # get phits-pmisses
    genes_in_pathway$phits_pmisses <- get_phits_pmisses(genes_in_pathway$phits, genes_in_pathway$pmisses)

    # append rows with NESrank
    genes_in_pathway<- merge(genes_in_pathway, pathways_significant, by.x="pathway_id", by.y="Pathway")
    genes_in_pathway
  }

  rugplots_data <- rugplots_data %>% select(pathway_id, NESrank, gene_id, rank, phits_pmisses)
  rugplots_data <- merge(rugplots_data, select(pathways_unique, Pathway, pvalue, FDR, qvalue), by.x="pathway_id", by.y="Pathway") %>% arrange(NESrank)
  rugplots_data
}
