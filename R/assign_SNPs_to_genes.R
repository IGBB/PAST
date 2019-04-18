# These aren't really global, but R CMD check
# doesn't understand how to deal with variable
# created by foreach loops
globalVariables("snp_chunk")
globalVariables("block")

#' Determine Linkage
#'
#' @param chunk A chunk of data to be processed
#' @param r_squared_cutoff The R^2 value to check against
#' @importFrom rlang .data
#' @import dplyr
#' @return Either the first unlinked SNP or a set of linked SNPs
determine_linkage <- function(chunk, r_squared_cutoff) {
  return_chunk = NULL
  for (i in seq_along(chunk)){
    data <- chunk[[i]]
    data <- data %>% dplyr::arrange(.data$Dist_bp)
    linked <- data %>% dplyr::filter(.data$R.2 >= r_squared_cutoff)
    unlinked <- data %>% dplyr::filter(.data$R.2 < r_squared_cutoff)

    if (nrow(unlinked) == nrow(data)) {
      return_chunk <- rbind(return_chunk, unlinked[1, ])
    } else {
      return_chunk <- rbind(return_chunk,
                            linked %>% filter(.data$R.2 >= r_squared_cutoff))
    }
  }
  return_chunk
}

#' Find representative SNP for a chunk of SNPs
#'
#' @param chunk A chunk of data to parse
#' @param r_squared_cutoff The R^2 value to check against when counting SNPs
#' @importFrom rlang .data
#' @import dplyr
#' @return A single SNP representing the whole chunk
find_representative_SNP <- function(chunk, r_squared_cutoff) {
  chunk <- chunk %>%
    dplyr::mutate(linked_snp_count = nrow(chunk),
                  linked_snp_count = ifelse(.data$R.2 < r_squared_cutoff,
                                            0,
                                            .data$linked_snp_count))
  if (chunk[1, ]$linked_snp_count > 1) {

    # count the number of negative/positive effects in each chunk
    negative <- sum(chunk$SNP2_effect < 0)
    positive <- sum(chunk$SNP2_effect > 0)

    # Find SNP with largest negative or positive effect
    if (positive > negative) {
      # sort in descending order
      chunk <- chunk %>%
        dplyr::arrange(desc(.data$SNP2_effect), .data$Dist_bp)

      # get top row
      chunk <- chunk[1, ]

    } else if (negative > positive) {

      # sort in ascending order
      chunk <- chunk %>%
        dplyr::arrange(.data$SNP2_effect, .data$Dist_bp)

      # get top row
      chunk <- chunk[1, ]

    } else if (negative == positive) {
      if (chunk$SNP1_effect[[1]] > 0) {

        # sort in descending order
        chunk <- chunk %>%
          dplyr::arrange(desc(.data$SNP2_effect))

        # get top row
        chunk <- chunk[1, ] %>%
          dplyr::mutate(linked_snp_count = -1)
      } else{
        # sort in ascending order
        chunk <- chunk %>%
          dplyr::arrange(.data$SNP2_effect)

        # get top row
        chunk <- chunk[1, ] %>%
          dplyr::mutate(linked_snp_count = -1)
      }
    }
  }
  chunk
}

#' Assign SNPs in a chunk to genes
#'
#' @param gff The GFF data for the chromosome being parsed
#' @param chunk The dataframe containing SNP data
#' @param window The search window around the SNPs
#' @importFrom rlang .data
#' @import dplyr
#' @return tagSNPs labeled with gene names
assign_chunk <- function(gff, chunk, window) {

  # set up data.frame of conditions to select tagSNP
  conditions <- chunk %>%
    dplyr::mutate(unlinked = .data$linked_snp_count == 0,
           single_linked = .data$linked_snp_count == 1,
           many_linked = .data$linked_snp_count > 1,
           problem_linked = .data$linked_snp_count == -1,
           effects_same_sign = .data$SNP1_effect * .data$SNP2_effect > 0,
           largest_effect = ifelse(abs(.data$SNP1_effect) >
                                     abs(.data$SNP2_effect),
                                   "SNP1",
                                   "SNP2"),
           effects_same_sign_use = ifelse(.data$SNP1_effect == .data$SNP2_effect,
                                          "SNP2",
                                          largest_effect),
           pvalue_equal = ifelse(SNP1_pvalue == SNP2_pvalue,
                                 TRUE,
                                 FALSE),
           pvalue_not_equal_use = ifelse(SNP1_pvalue > SNP2_pvalue,
                                         "SNP2",
                                         "SNP1"),
           single_pvalue_use = ifelse(pvalue_equal == FALSE,
                                      pvalue_not_equal_use,
                                      "PROBLEM"),
           single_use = ifelse(effects_same_sign == TRUE,
                               effects_same_sign_use,
                               pvalue_not_equal_use),
           distance_use = ifelse(.data$Dist_bp < window,
                                 effects_same_sign_use,
                                 "SNP2"),
           many_use = ifelse(effects_same_sign == TRUE,
                             distance_use,
                             "SNP2"),
           problem_use = ifelse(effects_same_sign == TRUE,
                                effects_same_sign_use,
                                "PROBLEM"),
           chromosome = .data$Locus,
           position = NA,
           effect = NA,
           p.value = NA,
           distance = .data$Dist_bp) %>%
    dplyr::select(.data$unlinked,
           .data$single_linked,
           .data$many_linked,
           .data$problem_linked,
           .data$single_use,
           .data$many_use,
           .data$problem_use,
           .data$Position1,
           .data$Position2,
           .data$SNP1_pvalue,
           .data$SNP1_effect,
           .data$SNP2_pvalue,
           .data$SNP2_effect,
           .data$chromosome,
           .data$position,
           .data$effect,
           .data$p.value,
           .data$distance,
           .data$linked_snp_count
    )

  # handle unlinked SNPs
  conditions <- rbind(conditions %>% dplyr::filter(.data$unlinked == FALSE),
                      conditions %>% dplyr::filter(.data$unlinked == TRUE) %>%
                        dplyr::mutate(position = .data$Position1,
                               effect = .data$SNP1_effect,
                               p.value = .data$SNP1_pvalue)
  ) %>%
    dplyr::select(-unlinked)

  # handle single linked SNPs
  conditions <- rbind(conditions %>%
                        dplyr::filter(.data$single_linked == FALSE),
                      conditions %>%
                        dplyr::filter(.data$single_linked == TRUE) %>%
                        dplyr::mutate(position = ifelse(.data$single_use == "SNP1",
                                                 .data$Position1,
                                                 .data$Position2),
                               effect = ifelse(.data$single_use == "SNP1",
                                               .data$SNP1_effect,
                                               .data$SNP2_effect),
                               p.value = ifelse(.data$single_use == "SNP1",
                                                .data$SNP1_pvalue,
                                                .data$SNP2_pvalue))
  ) %>%
    dplyr::select(-single_linked,
           -single_use)

  # handle multiply-linked SNPs
  conditions <- rbind(conditions %>%
                        dplyr::filter(.data$many_linked == FALSE),
                      conditions %>%
                        dplyr::filter(.data$many_linked == TRUE) %>%
                        dplyr::mutate(position = ifelse(.data$many_use == "SNP1",
                                                 .data$Position1,
                                                 .data$Position2),
                               effect = ifelse(.data$many_use == "SNP1",
                                               .data$SNP1_effect,
                                               .data$SNP2_effect),
                               p.value = ifelse(.data$many_use == "SNP1",
                                                .data$SNP1_pvalue,
                                                .data$SNP2_pvalue))
  ) %>%
    dplyr::select(-many_linked,
           -many_use)

  # handle problem SNPs
  tagSNPs <- rbind(conditions %>%
                     dplyr::filter(.data$problem_linked == FALSE),
                   conditions %>%
                     dplyr::filter(.data$problem_linked == TRUE) %>%
                     dplyr::mutate(position = ifelse(.data$problem_use == "SNP1",
                                              .data$Position1,
                                              .data$Position2),
                            effect = ifelse(.data$problem_use == "SNP1",
                                            .data$SNP1_effect,
                                            .data$SNP2_effect),
                            p.value = ifelse(.data$problem_use == "SNP1",
                                             .data$SNP1_pvalue,
                                             .data$SNP2_pvalue))
  ) %>%
    dplyr::select(.data$chromosome,
           .data$position,
           .data$effect,
           .data$p.value,
           .data$distance,
           .data$linked_snp_count) %>%
    dplyr::mutate(window_start = .data$position - window,
           window_end = .data$position + window) %>%
    dplyr::arrange(position)

  # assign genes
  inner_join(tagSNPs, gff, by = c("chromosome" = "seqid")) %>%
    dplyr::filter(
      (
        .data$window_start >= .data$start &
          .data$window_end <= .data$end
      ) |
        (
          .data$start >= .data$window_start & .data$end <= .data$window_end
        ) |
        (
          .data$window_start >= .data$start &
            .data$window_start <= .data$end
        ) |
        (
          .data$window_end >= .data$start &
            .data$window_end <= .data$end
        )
    ) %>%
    dplyr::select(-c(
      .data$distance,
      .data$window_start,
      .data$window_end,
      .data$source,
      .data$type,
      .data$start,
      .data$end,
      .data$score,
      .data$strand,
      .data$phase,
      .data$ID,
      .data$biotype
    )) %>%
    dplyr::mutate(name = .data$Name,
           Name = NULL)
}

#' Find the SNP-gene assignment that represents SNPs assigned to a gene
#'
#' @param chunk A chunk of gene assignments
#' @importFrom rlang .data
#' @import dplyr
#' @return A single SNP-gene assignment representing all SNPS assigned to the
#'   same gene
#' to a gene
find_representative_SNP_gene_pairing <- function(chunk) {
  neg_genes <- chunk %>%
    dplyr::filter(.data$effect < 0) %>%
    dplyr::arrange(.data$effect, .data$p.value) %>%
    dplyr::mutate(linked_snp_count = ifelse(.data$linked_snp_count <= 0,
                                            1,
                                            .data$linked_snp_count))
  pos_genes <- chunk %>%
    dplyr::filter(.data$effect > 0) %>%
    dplyr::arrange(desc(.data$effect), .data$p.value) %>%
    dplyr::mutate(linked_snp_count = ifelse(.data$linked_snp_count <= 0,
                                            1,
                                            .data$linked_snp_count))
  negative <- sum(neg_genes$linked_snp_count)
  positive <- sum(pos_genes$linked_snp_count)

  if (positive > negative) {
    chunk <- chunk %>%
      dplyr::arrange(desc(.data$effect), desc(.data$p.value))
    representative <- chunk[1, ]
  } else if (negative > positive) {
    chunk <- chunk %>%
      dplyr::arrange(.data$effect, .data$p.value)
    representative <- chunk[1, ]
  } else if (positive == negative) {
    pos_max <- pos_genes[1, ]
    neg_max <- neg_genes[1, ]

    if (pos_max$effect > abs(neg_max$effect)) {
      representative <- pos_max
    }
    else{
      representative <- neg_max
    }
  }
  representative %>% dplyr::mutate(linked_snp_count = negative + positive)
}

#' Assign SNPs to genes
#'
#' @param gwas_data Merged association and effects data from merge_data()
#' @param LD Linkage disequilibrium data from parse_LD()
#' @param gff_file The path to a GFF file
#' @param window The search window for genes around the SNP
#' @param r_squared_cutoff The R^2 value used to determine SNP significance
#' @param num_cores The number of cores to use in parallelizing PAST
#' @importFrom rlang .data
#' @importFrom rtracklayer readGFF
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @return A dataframe of genes from the SNP data
#' @export
#'
#' @examples
#' demo_association_file = system.file("extdata", "association.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_GWAS_data(demo_association_file, demo_effects_file)
#' demo_linkage_disequilibrium_file = system.file("extdata","LD.txt.xz",
#' package = "PAST", mustWork = TRUE)
#' LD <- load_LD(demo_linkage_disequilibrium_file)
#' demo_genes_file = system.file("extdata", "genes.gff",
#'   package = "PAST", mustWork = TRUE)
#' genes <-assign_SNPs_to_genes(gwas_data, LD, demo_genes_file, 1000, 0.8, 2)
assign_SNPs_to_genes <-
  function(gwas_data,
           LD,
           gff_file,
           window,
           r_squared_cutoff,
           num_cores) {

    full_gff <- rtracklayer::readGFF(gff_file) %>% filter(.data$type == "gene")
    chromosomes <- as.character.factor(
      full_gff %>%
        dplyr::select(.data$seqid) %>%
        dplyr::arrange(.data$seqid) %>%
        unique() %>%
        dplyr::pull(.data$seqid)
    )

    cl <- parallel::makeCluster(num_cores)
    registerDoParallel(cl)

    all_genes <- NULL

    # UP/DOWNSTREAM LOOP
    for (i in seq_along(c(1, 2))) {

      # BEGIN PROCESSING BY CHROMOSOMES LOOP
      for (chromosome in names(LD)) {

        if (chromosome %in% chromosomes) {
          temp_data <-
            LD[[chromosome]] %>%
            dplyr::mutate(Marker1 = paste0("S",
                                           .data$Locus,
                                           "_",
                                           .data$Position1),
                          Marker2 = paste0("S",
                                           .data$Locus,
                                           "_",
                                           .data$Position2)) %>%
            dplyr::arrange(.data$Position1)

          if (i == 2) {
            temp_data <-
              temp_data %>%
              dplyr::mutate(
                temp_p2 = .data$Position1,
                temp_s2 = .data$Site1,
                temp_Marker2 = .data$Marker1,
                Position1 = .data$Position2,
                Site1 = .data$Site2,
                Marker1 = .data$Marker2
              ) %>%
              dplyr::mutate(
                Site2 = .data$temp_s2,
                Position2 = .data$temp_p2,
                Marker2 = .data$temp_Marker2,
                temp_s2 = NULL,
                temp_p2 = NULL,
                temp_Marker2 = NULL
              ) %>%
              dplyr::arrange(.data$Site1)
          }
          temp_data_list <- split(temp_data, temp_data$Position1)
          temp_data_list_split <- split(temp_data_list,
                                        ceiling(seq_along(temp_data_list)
                                                /num_cores))

          temp_data <-
            foreach(
              data = temp_data_list_split,
              .combine = rbind,
              .packages = "dplyr"
            ) %dopar% {
              determine_linkage(data, r_squared_cutoff)
            }

          gwas_data_for_chromosome <-
            dplyr::filter(gwas_data, .data$Chr == as.integer(chromosome))

          # look up p-value and effect data for SNP1
          temp_data <-
            merge(temp_data, gwas_data_for_chromosome,
                  by.x = "Marker1",
                  by.y = "Marker") %>%
            dplyr::mutate(SNP1_pvalue = .data$p,
                          SNP1_effect = .data$Effect.x) %>%
            dplyr::select(
              .data$Locus,
              .data$Position1,
              .data$Position2,
              .data$Site1,
              .data$Site2,
              .data$Dist_bp,
              .data$R.2,
              .data$Marker1,
              .data$SNP1_pvalue,
              .data$SNP1_effect,
              .data$Marker2
            )

          # look up p-value and effect data for SNP2
          temp_data <-
            merge(temp_data, gwas_data_for_chromosome,
                  by.x = "Marker2",
                  by.y = "Marker") %>%
            dplyr::mutate(SNP2_pvalue = .data$p,
                          SNP2_effect = .data$Effect.x) %>%
            dplyr::select(
              .data$Locus,
              .data$Position1,
              .data$Position2,
              .data$Site1,
              .data$Site2,
              .data$Dist_bp,
              .data$R.2,
              .data$Marker1,
              .data$SNP1_pvalue,
              .data$SNP1_effect,
              .data$Marker2,
              .data$SNP2_pvalue,
              .data$SNP2_effect
            ) %>%
            dplyr::arrange(.data$Position1)

          index <- c(0, cumsum(abs(diff(temp_data$Site2)) > 1))
          temp_data_list <- split(temp_data, paste(temp_data$Position1, index))

          temp_data <-
            foreach(
              data = temp_data_list,
              .combine = rbind,
              .packages = "dplyr"
            ) %dopar% {
              find_representative_SNP(data, r_squared_cutoff)
            }

          split <- 1000
          snp_list <-
            split(temp_data,
                  rep(seq_along(split),
                      length.out = nrow(temp_data),
                      each = ceiling(nrow(temp_data) / split)
                  ))

          # subset gff to only handle this chromosome
          gff <- dplyr::filter(full_gff, .data$seqid == chromosome) %>%
            mutate(seqid = as.character.factor(.data$seqid))

          # get genes in parallel
          chromosome_genes <-
            foreach(
              snp_chunk = snp_list,
              .combine = rbind,
              .packages = "dplyr"
            ) %dopar% {
              assign_chunk(gff, snp_chunk, window)
            }

          all_genes <- rbind(all_genes, chromosome_genes)
        }
      }
    }

    group_by_gene <- split(all_genes, f = all_genes$name)
    representative_genes <-
      foreach(
        chunk = group_by_gene,
        .combine = rbind,
        .packages = c("dplyr", "PAST")
      ) %dopar% {
        find_representative_SNP_gene_pairing(chunk)
      }

    stopCluster(cl)

    representative_genes
  }
