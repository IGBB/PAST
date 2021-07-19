#' Find Pathway Significance
#'
#' @param genes genes from assign_SNPs_to_genes()
#' @param pathways_file a file containing the pathway IDs, their names, and the
#'   genes in the pathway with one row per gene for every pathway
#' @param gene_number_cutoff a cut-off for the minimum number of genes in a
#'   pathway; pathways that don't satisfy this requirement are removed from analysis
#' @param analysis_mode increasing/decreasing, where increasing analysis mode 
#' tests whether pathways are associated with an increase in the measured trait 
#' and decreasing analysis mode tests whether pathways are associated with a 
#' decrease in the measured trait 
#' @param permutations How many times to sample the effects data during random
#'   sampling
#' @import data.table
#' @import qvalue
#' @return enrichment data that indicates which pathways are significantly
#' associated with the trait (via p-value/q-value) and contains information for
#' plotting running enrichment score for genes in a pathway
#' @export
#'
#' @examples
#' example("assign_SNPs_to_genes")
#' pathways = system.file("extdata", "pathways.txt.gz", package = "PAST", 
#' mustWork = TRUE)
#' enrichment_data <- find_pathway_significance(genes, pathways, 5, "increasing", 1000)
find_pathway_significance <-
  function(genes,
           pathways_file,
           gene_number_cutoff = 5,
           analysis_mode,
           permutations = 1000) {
    
    message(data.table::getDTthreads(verbose = getOption("datatable.verbose")))
    
    # Read pathways file.
    # PAST technically accepts a data.table of pathways or a path to a file
    #   containing pathways information.
    # However, PAST only accepts a data.table because doing so makes working
    #   with PAST Shiny easier. Users should pass a path.
    if (data.table::is.data.table(pathways_file) == FALSE) {
      pathways <- data.table::fread(pathways_file, quote = "")  
    } else {
      pathways <- data.table::copy(pathways_file)
    }
    
    # Order the pathways data and drop duplicate entries.
    setkey(pathways)
    pathways <- unique(pathways)
    
    # PAST computes pathway significance based on a null distribution generated
    #   from permutations of the true effects and a user-provided number of
    #   permutations.
    # Effects are sorted and ranked based on the provided analysis mode.
    # - "increasing" sorts effects from most positive to most negative.
    # - "decreasing" sorts effects from most negative to most positive and 
    #      converts the effects to the absolute value of the effects.
    effects = data.table::copy(genes[, .(name, effect)])
    data.table::setnames(effects, "name", "gene")
    data.table::setnames(effects, "effect", "effect0")
    if (analysis_mode == "increasing") {
      setorder(effects, -effect0)
    } else if (analysis_mode == "decreasing") {
      setorder(effects, effect0)
      effects[, paste("effect", column, sep="") := abs(paste("effect", column, sep=""))]
    } else {
      stop("Analysis mode must be \"increasing\" or \"decreasing\".")
    }
    effects[, rank0 := as.numeric(.I)]
    columns = as.character(1:permutations)
    for (column in columns) {
      data.table::set(effects, j = paste("effect", column, sep=""), value = sample(genes[,effect]))
      if (analysis_mode == "increasing") {
        setorderv(effects, paste("effect", column, sep=""), -1)
      } else if (analysis_mode == "decreasing") {
        setorderv(effects, paste("effect", column, sep=""))
        effects[, paste("effect", column, sep="") := abs(paste("effect", column, sep=""))]
      } else {
        stop("Analysis mode must be \"increasing\" or \"decreasing\".")
      }
      effects[, paste("rank", column, sep = "") := as.numeric(.I)]
    }
    
    # Find enrichment score for permutations.
    # Count the number of genes overall.
    # Filter pathways such that only genes with an effect are kept.
    number_of_genes = genes[,.N]
    enrichment_scores = copy(pathways[ 
      effects[gene %in% pathways[,gene_id]],
      on = .(gene_id = gene)
    ])
    
    # Find enrichment score for each permutation of the effects.
    # Show progress by pathway_id.
    # Method described in DOI: 10.1073/pnas.0506580102
    number_of_pathways = data.table::uniqueN(enrichment_scores$pathway_id)
    progress <- utils::txtProgressBar(min = 0, max = number_of_pathways, style = 3)
    enrichment_scores[, (paste0("enrichment_score", c(0:permutations))) := lapply(c(0:permutations), function(column) {
      setorderv(.SD, paste0("rank", column))
      
      # Calculate factors.
      # Factors are calculated by subtracting 1 + the rank of the previous gene
      #   from the rank of the current gene.
      # The cumulative sum of the factors per row is calculated and added to the
      #   previously calculated factor.
      factor = unlist(.SD[,paste0("rank", column), with = F] - data.table::shift(.SD[,paste0("rank", column), with = F], fill = 0) - 1)
      factor = factor + data.table::shift(cumsum(factor), fill = 0)
      
      # The factor is divided by the number of genes that are not in the pathway
      #   to calculate the pmiss.
      pmiss = unlist(factor / (number_of_genes - .N))
      
      # The effect of each gene in the pathway is divided by the sum of the
      #   absolute value of every effect in the pathway.
      # A column is created to track the phit_base of the previous gene in the
      #   pathway.
      # The cumulative sum of the lagged phits per row is calculated and added 
      # to the phit.
      phit = unlist(abs(.SD[, paste0("effect", column), with = F]) / sum(abs(.SD[, paste0("effect", column), with = F])))
      phit = phit + cumsum(data.table::shift(phit, fill = 0))
      
      # The progress bar is updated.
      # The pmiss is subtracted from the phit to calculate the enrichment score
      #   per gene.
      # The maximum enrichment score from the enrichment scores of all genes in
      #   the pathway is returned.
      utils::setTxtProgressBar(progress, .GRP)
      max(phit - pmiss)
    }),
    by = pathway_id,
    .SDcols = names(enrichment_scores)[3:length(names(enrichment_scores))]]
    
    # enrichment_scores_safe = copy(enrichment_scores)
    # enrichment_scores_safe = readRDS("enrichment_scores.Rds")
    # saveRDS(enrichment_scores_safe, "enrichment_scores.Rds")
    # enrichment_scores = copy(pathways_safe)
    
    # Drop pathways without enough genes to meet the cutoff filter.
    # Drop rank, effect, and genes.
    # Take the first row by pathway.
    # Rename enrichment_score0 to enrichment_score_observed.
    enrichment_scores = enrichment_scores[, if(.N > gene_number_cutoff) .SD, by=pathway_id]
    enrichment_scores[, c("gene_id", paste0("rank", c(0:permutations)), paste0("effect", c(0:permutations))) := NULL]
    enrichment_scores = enrichment_scores[, .SD[1L], by = pathway_id]
    setnames(enrichment_scores, "enrichment_score0", "enrichment_score_observed")
    enrichment_scores = data.table::melt(enrichment_scores, 
                                         id.vars = c("pathway_id", 
                                                     "pathway_name", 
                                                     "enrichment_score_observed"),
                                         measure.vars = paste0("enrichment_score", c(1:permutations)))
    setkey(enrichment_scores, pathway_id)
    enrichment_scores[,permutation_mean := mean(value)]
    enrichment_scores[,permutation_sd := sd(value)]
    enrichment_scores = dcast(enrichment_scores, pathway_id + pathway_name + enrichment_score_observed + permutation_mean + permutation_sd ~ variable)
    enrichment_scores[, paste0("enrichment_score", c(1:permutations)) := NULL]
    enrichment_scores[, normalized_enrichment_score := (enrichment_score_observed - permutation_mean) / permutation_sd]
    enrichment_scores[, `p-value` := 1-pnorm(normalized_enrichment_score)]
    enrichment_scores[, `q-value` := qvalue::qvalue(`p-value`, fdr.level=0.05, lambda = 0)$qvalues]
    enrichment_scores[, c("permutation_mean", "permutation_sd", "normalized_enrichment_score") := NULL]
    
    # Calculate the running enrichment scores per gene using the same method as
    #   used for calculating the enrichment score for for a pathway.
    # Keep per-gene enrichment scores, rather than returning the max.
    running_enrichment_scores = copy(pathways[ 
      effects[gene %in% pathways[,gene_id], .(gene, effect0, rank0)],
      on = .(gene_id = gene)
    ])
    setnames(running_enrichment_scores, "effect0", "effect")
    setnames(running_enrichment_scores, "rank0", "rank")
    setorder(running_enrichment_scores, rank)
    running_enrichment_scores[, enrichment_score := {
      
      # Calculate factors.
      # Factors are calculated by subtracting 1 + the rank of the previous gene
      #   from the rank of the current gene.
      # The cumulative sum of the factors per row is calculated and added to the
      #   previously calculated factor.
      factor = unlist(.SD[, rank] - data.table::shift(.SD[, rank], fill = 0) - 1)
      factor = factor + data.table::shift(cumsum(factor), fill = 0)
      
      # The factor is divided by the number of genes that are not in the pathway
      #   to calculate the pmiss.
      pmiss = unlist(factor / (number_of_genes - .N))
      
      # The effect of each gene in the pathway is divided by the sum of the
      #   absolute value of every effect in the pathway.
      # A column is created to track the phit_base of the previous gene in the
      #   pathway.
      # The cumulative sum of the lagged phits per row is calculated and added 
      # to the phit.
      phit = unlist(abs(.SD[, effect]) / sum(abs(.SD[, effect])))
      phit = phit + cumsum(data.table::shift(phit, fill = 0))
      
      # The pmiss is subtracted from the phit to calculate the enrichment score
      #   per gene.
      phit - pmiss
    },
    by = pathway_id,
    .SDcols = names(running_enrichment_scores)]
    
    # Drop pathways without enough genes to meet the cutoff filter.
    # Combine the enrichment score data for pathways with the running enrichment
    #   score data for genes.
    running_enrichment_scores = running_enrichment_scores[, if(.N > gene_number_cutoff) .SD, by=pathway_id]
    enrichment_data = enrichment_scores[
      running_enrichment_scores[, .(pathway_id, gene_id, rank, enrichment_score)],
      on = .(pathway_id)]
    data.table::setorder(enrichment_data, `q-value`, pathway_id)
    return(enrichment_data)
  }