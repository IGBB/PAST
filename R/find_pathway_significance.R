get_factors <- function(gene_ranks) {
  factors <- vector("logical", length(gene_ranks))
  factors[1] <- gene_ranks[1] - 1
  for (i in 2:length(gene_ranks)) {
    factors[i] <- gene_ranks[i] - gene_ranks[i - 1] - 1 + factors[i - 1]
  }
  factors
}

get_pmisses <- function(temp_data, genes_in_pathway, factors) {
  pmisses <- factors / (nrow(temp_data) - nrow(genes_in_pathway))
  pmisses
}

get_phits <- function(gene_effects) {
  NR <- sum(abs(gene_effects))
  absolute <- abs(gene_effects)
  phits <- vector("logical", length(gene_effects))
  phits[1] <- (absolute[1] / NR)
  for (i in 2:length(gene_effects)) {
    phits[i] <- absolute[i] / NR + phits[i - 1]
  }
  phits
}

get_running_enrichment_score <- function(phits, pmisses) {
  running_enrichment_score <- phits - pmisses
  running_enrichment_score
}

find_max <- function(running_enrichment_score) {
  running_enrichment_score <- sort(running_enrichment_score, decreasing = TRUE)
  max <- running_enrichment_score[[1]]
  max
}

# This isn't really global, but R CMD check
# doesn't understand how to deal with variable
# created by foreach loops
globalVariables("pathway")

#' Find Pathway Significance
#'
#' @param genes Genes from assign_SNPs_to_genes()
#' @param pathways_file A file containing the pathway IDs, their names, and the
#'   genes in the pathway
#' @param gene_number_cutoff A cut-off for the minimum number of genes in a
#'   pathway
#' @param mode increasing/decreasing
#' @param sample_size How many times to sample the effects data during random
#'   sampling
#' @param num_cores The number of cores to use in parallelizing PAST
#' @importFrom rlang .data
#' @importFrom stats pnorm
#' @importFrom stats sd
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @import iterators
#' @import parallel
#' @import qvalue
#' @return Rugplots data
#' @export
#'
#' @examples
#' example("assign_SNPs_to_genes")
#' demo_pathways_file = system.file("extdata", "pathways.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' rugplots_data <- find_pathway_significance(genes, demo_pathways_file, 5,
#'   "increasing", 1000, 2)
find_pathway_significance <-
  function(genes,
           pathways_file,
           gene_number_cutoff = 5,
           mode,
           sample_size = 1000,
           num_cores) {
    # load pathways
    pathways <-
      read.table(pathways_file,
                 sep = "\t",
                 header = TRUE,
                 quote = "")

    # sample to create 1000 random distributions
    effects <- genes %>% dplyr::select(.data$name, .data$effect)
    effects <-
      cbind(effects, vapply(seq_len(sample_size),
                            function(i) sample(effects$effect),
                            FUN.VALUE = double(nrow(effects))))

    pathways_unique <- unique(select(pathways, .data$pathway_id))
    pathways_unique[] <- lapply(pathways_unique, as.character)

    # process sample columns
    cl <- parallel::makeCluster(num_cores, outfile = "")
    registerDoParallel(cl)
    column_observations <-
      foreach(
        i = iter(2:(sample_size + 2)),
        .combine = cbind,
        .packages = c("dplyr", "PAST", "foreach", "iterators")
      ) %dopar% {
        temp_data <-
          data.frame(matrix("NA", ncol = 2, nrow = nrow(effects)))
        colnames(temp_data) <- c("gene", "effect")
        temp_data$effect <- effects[, i]
        temp_data$gene <- effects[, 1]

        if (mode == "decreasing") {
          temp_data <- temp_data %>%
            dplyr::arrange(.data$effect) %>%
            dplyr::mutate(rank = row_number(),
                          effect = abs(.data$effect))
        } else if (mode == "increasing") {
          temp_data <- temp_data %>%
            dplyr::arrange(desc(.data$effect)) %>%
            dplyr::mutate(rank = row_number())
        } else {
          stop("Incorrect mode.")
        }

        column_observation <- data.frame()

        foreach(pathway = iter(pathways_unique$pathway_id, by = "row")) %do% {
          genes_in_pathway <-
            dplyr::filter(pathways, pathways$pathway_id == pathway) %>%
            filter(!(is.na(gene_id)))

          ## get ranks and effects and sort by rank
          genes_in_pathway <-
            merge(genes_in_pathway,
                  temp_data,
                  by.x = "gene_id",
                  by.y = "gene") %>%
            dplyr::arrange(.data$rank) %>% unique()

          # check cutoff
          if (nrow(genes_in_pathway) >= gene_number_cutoff) {

            # get factors using rank
            factors <- get_factors(genes_in_pathway$rank)

            # get pmisses
            pmisses <- get_pmisses(temp_data, genes_in_pathway, factors)

            # get phits
            phits <- get_phits(genes_in_pathway$effect)

            # get phits-pmisses
            running_enrichment_score <- get_running_enrichment_score(phits, pmisses)
            find_max(running_enrichment_score)
            # store max phit_pmisses
            column_observation <-
              rbind(column_observation, find_max(running_enrichment_score))
          } else {
            column_observation <- rbind(column_observation, NA)
          }
        }
        column_observation
      }

    stopCluster(cl)

    pathways_unique <- cbind(pathways_unique, column_observations)
    colnames(pathways_unique) <-
      c("Pathway", "ES_Observed", seq_len(sample_size))
    colnames(pathways_unique)[3:(sample_size + 2)] <-
      paste0("ES", colnames(pathways_unique)[3:(sample_size + 2)])
    pathways_unique <- pathways_unique %>%
      dplyr::filter(!is.na(.data$ES_Observed))
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(permutation_mean =
                      apply(pathways_unique[, 3:(sample_size + 2)], 1, mean))
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(permutation_standard_deviation =
                      apply(pathways_unique[, 3:(sample_size + 2)], 1, sd))
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(
        NES_Observed = (.data$ES_Observed - .data$permutation_mean) / .data$permutation_standard_deviation
      )
    for (i in seq_along(sample_size)) {
      pathways_unique[i + 2] <-
        (pathways_unique[i + 2] - pathways_unique[sample_size + 3]) / pathways_unique[sample_size + 4]
    }
    pathways_unique <- pathways_unique %>%
      dplyr::mutate(pvalue = 1-pnorm(.data$NES_Observed))
    pathways_unique <-dplyr::arrange(pathways_unique, .data$pvalue)
    pathways_unique <- pathways_unique %>%
      dplyr::select(.data$Pathway,
                    .data$ES_Observed,
                    .data$NES_Observed,
                    .data$pvalue) %>%
      dplyr::mutate(qvalue = qvalue(
        pathways_unique$pvalue,
        lambda = 0,
        fdr.level = 0.05
      )$qvalues)

    pathways_significant <- pathways_unique %>%
      dplyr::select(.data$Pathway,
                    .data$NES_Observed,
                    .data$pvalue,
                    .data$qvalue) %>%
      dplyr::mutate(pathway_number = row_number(), NESrank = NULL)

    temp_data <-
      data.frame(matrix("NA", ncol = 2, nrow = nrow(effects)))
    colnames(temp_data) <- c("Gene", "Effect")
    temp_data$Effect <- effects[, 2]
    temp_data$Gene <- effects[, 1]
    colnames(temp_data) <- c("Gene", "Effect")
    if (mode == "decreasing") {
      temp_data <- temp_data %>%
        dplyr::arrange(.data$Effect) %>%
        dplyr::mutate(rank = row_number(),
                      effect = abs(.data$Effect))
    } else if (mode == "increasing") {
      temp_data <- temp_data %>%
        dplyr::arrange(desc(.data$Effect)) %>%
        dplyr::mutate(rank = row_number())
    } else {
      stop("Incorrect mode.")
    }

    rugplots_data <- NULL
    rugplots_data <-
      foreach(
        pathway = iter(pathways_significant$Pathway, by = "row"),
        .combine = rbind
      ) %do% {
        genes_in_pathway <-
          dplyr::filter(pathways, pathways$pathway_id == pathway)

        ## get ranks and effects and sort by rank
        genes_in_pathway <-
          merge(genes_in_pathway,
                temp_data,
                by.x = "gene_id",
                by.y = "Gene") %>%
          dplyr::arrange(rank)

        # get factors using rank
        genes_in_pathway$factors <-
          get_factors(genes_in_pathway$rank)

        # get pmisses
        genes_in_pathway$pmisses <-
          get_pmisses(temp_data, genes_in_pathway, genes_in_pathway$factors)

        # get phits
        genes_in_pathway$phits <- get_phits(genes_in_pathway$Effect)

        # get phits-pmisses
        genes_in_pathway$running_enrichment_score <-
          get_running_enrichment_score(genes_in_pathway$phits, genes_in_pathway$pmisses)

        # append rows with NESrank
        genes_in_pathway <-
          merge(genes_in_pathway,
                pathways_significant,
                by.x = "pathway_id",
                by.y = "Pathway")
        genes_in_pathway
      }
    rugplots_data <- merge(
      rugplots_data %>%
        dplyr::select(
          .data$pathway_id,
          .data$pathway_number,
          .data$gene_id,
          .data$rank,
          .data$running_enrichment_score
        ),
      pathways %>%
        dplyr::select(.data$pathway_id, .data$pathway_name) %>%
        unique(),
      by = "pathway_id"
    )
    rugplots_data <- merge(
      rugplots_data,
      dplyr::select(
        pathways_unique,
        .data$Pathway,
        .data$pvalue,
        .data$qvalue
      ),
      by.x = "pathway_id",
      by.y = "Pathway"
    ) %>%
      dplyr::arrange(.data$pathway_number)
    rugplots_data
  }