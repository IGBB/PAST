#' Read a GFF file
#'
#' @param gff_file The path to the GFF file
#' @importFrom rlang .data
#' @import dplyr
#' @return A dataframe containing the relevant GFF data - chromosome, gene name, start and stop locations
read_gff <- function(gff_file) {
  gff <- read.table(
    gff_file,
    comment.char = "#",
    col.names =
      c(
        "chr",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
      )
  )  %>%
    dplyr::select(.data$chr, .data$attributes, .data$start, .data$end) %>%
    dplyr::mutate(tempcol = stringr::str_split(attributes, ";Name=")) %>%
    rowwise() %>%
    dplyr::mutate(
      name = gsub("ID=", "", unlist(.data$tempcol)[1]),
      tempcol = NULL,
      attributes = NULL,
      chr = gsub("chr", "", .data$chr)
    )
  gff
}

#' Parse Chunks of Data for Significant SNPs
#'
#' @param data A chunk of data to be processed
#' @param r_squared_cutoff The R^2 value to check against
#' @importFrom rlang .data
#' @import dplyr
#' @return Either the first unlinked SNP or a set of linked SNPs
parse_temp <- function(data, r_squared_cutoff) {
  data <- data %>% dplyr::arrange(.data$Dist_bp)
  linked <- data %>% dplyr::filter(.data$R.2 >= r_squared_cutoff)
  unlinked <- data %>% dplyr::filter(.data$R.2 < r_squared_cutoff)

  if (nrow(unlinked) == nrow(data)) {
    return <- unlinked[1, ]
  } else {
    return <- linked %>% filter(.data$R.2 >= r_squared_cutoff)
  }
  return
}

#' Parse blocked SNPs into a single representative SNP
#'
#' @param block A block of data to parse
#' @param r_squared_cutoff The R^2 value to check against when counting SNPs
#' @importFrom rlang .data
#' @import dplyr
#' @return A single SNP representing the whole block
parse_chunk <- function(block, r_squared_cutoff) {
  block <- block %>% dplyr::mutate(linked_snp_count = nrow(block))
  if (block[1, ]$linked_snp_count > 1) {
    # count the number of negative/positive effects in each block
    negative <- sum(block$SNP2_effect < 0)
    positive <- sum(block$SNP2_effect > 0)

    # Find SNP with largest negative or positive effect
    if (positive > negative) {
      # sort in descending order
      block <-
        block %>% dplyr::arrange(desc(.data$SNP2_effect), .data$Dist_bp)

      # get top row
      block <- block[1, ]

    } else if (negative > positive) {
      # sort in ascending order
      block <-
        block %>% dplyr::arrange(.data$SNP2_effect, .data$Dist_bp)

      # get top row
      block <- block[1, ]

    } else if (negative == positive) {
      if (block$SNP1_effect[[1]] > 0) {
        # sort in descending order
        block <- block %>% dplyr::arrange(desc(.data$SNP2_effect))

        # get top row
        block <- block[1, ] %>% dplyr::mutate(linked_snp_count = -1)
      } else{
        # sort in ascending order
        block <- block %>% dplyr::arrange(.data$SNP2_effect)

        # get top row
        block <- block[1, ] %>% dplyr::mutate(linked_snp_count = -1)
      }
    }
  } else {
    block <- block %>%
      dplyr::mutate(linked_snp_count = ifelse(.data$R.2 < r_squared_cutoff, 0, 1))
  }
  block
}

#' Select the SNP with the largest effect
#' from all SNPs assigned to a gene
#'
#' @param untagged_genes A chunk of untagged genes
#' @importFrom rlang .data
#' @import dplyr
#' @return A single SNP representing all SNPS assigned
#' to a gene
select_gene_from_block <- function(untagged_genes) {
  neg_genes <- untagged_genes %>%
    dplyr::filter(.data$Effect < 0) %>%
    dplyr::arrange(.data$Effect, .data$P_value) %>%
    dplyr::mutate(linked_snp_count = ifelse(.data$linked_snp_count <= 0, 1, .data$linked_snp_count))
  pos_genes <- untagged_genes %>%
    dplyr::filter(.data$Effect > 0) %>%
    dplyr::arrange(desc(.data$Effect), .data$P_value) %>%
    dplyr::mutate(linked_snp_count = ifelse(.data$linked_snp_count <= 0, 1, .data$linked_snp_count))
  negative <- sum(neg_genes$linked_snp_count)
  positive <- sum(pos_genes$linked_snp_count)

  if (positive > negative) {
    untagged_genes <-
      untagged_genes %>% dplyr::arrange(desc(.data$Effect), desc(.data$P_value))
    tag_snp <- untagged_genes[1, ]
  } else if (negative > positive) {
    untagged_genes <-
      untagged_genes %>% dplyr::arrange(.data$Effect, .data$P_value)
    tag_snp <- untagged_genes[1, ]
  } else if (positive == negative) {
    pos_max <- pos_genes[1, ]
    neg_max <- neg_genes[1, ]

    if (pos_max$Effect > abs(neg_max$Effect)) {
      tag_snp <- pos_max
    }
    else{
      tag_snp <- neg_max
    }
  }
  tag_snp %>% dplyr::mutate(linked_snp_count = negative + positive)
}

# These aren't really global, but R CMD check
# doesn't understand how to deal with variable
# created by foreach loops
globalVariables("snp_chunk")
globalVariables("block")

#' Find genes within a window of a SNP
#'
#' @param gff The GFF data for the chromosome being parsed
#' @param snp_df The dataframe containing SNP data
#' @param window The search window around the SNPs
#' @importFrom rlang .data
#' @import dplyr
#' @return SNPs labeled with gene names
find_genes <- function(gff, snp_df, window) {
  no_linked <- snp_df %>%
    dplyr::filter(.data$linked_snp_count == 0) %>%
    dplyr::mutate(Position = .data$Position1)
  one_linked <- snp_df %>%
    dplyr::filter(.data$linked_snp_count == 1) %>%
    dplyr::mutate(Position = ifelse(
      .data$SNP1_effect < 0 & .data$SNP2_effect < 0,
      ifelse(
        .data$SNP1_effect != .data$SNP2_effect,
        ifelse(
          .data$SNP1_effect < .data$SNP2_effect,
          .data$Position1,
          .data$Position2
        ),
        .data$Position2
      ),
      ifelse(
        .data$SNP1_effect > 0 & .data$SNP2_effect > 0,
        ifelse(
          .data$SNP1_effect != .data$SNP2_effect,
          ifelse(
            .data$SNP1_effect > .data$SNP2_effect,
            .data$Position1,
            .data$Position2
          ),
          .data$Position2
        ),
        ifelse(
          .data$SNP1_pval != .data$SNP2_pval,
          ifelse(
            .data$SNP1_pval < .data$SNP2_pval,
            .data$Position1,
            .data$Position2
          ),
          "problem"
        )
      )
    ))

  more_linked <- snp_df %>%
    dplyr::filter(.data$linked_snp_count > 1) %>%
    dplyr::mutate(Position = ifelse(
      .data$SNP1_effect < 0 & .data$SNP2_effect < 0,
      ifelse(
        .data$Dist_bp < window,
        ifelse(
          .data$SNP1_effect != .data$SNP2_effect,
          ifelse(
            .data$SNP1_effect < .data$SNP2_effect,
            .data$Position1,
            .data$Position2
          ),
          .data$Position2
        ),
        .data$Position2
      ),
      ifelse(
        .data$SNP1_effect > 0 & .data$SNP2_effect > 0,
        ifelse(
          .data$Dist_bp < window,
          ifelse(
            .data$SNP1_effect != .data$SNP2_effect,
            ifelse(
              .data$SNP1_effect > .data$SNP2_effect,
              .data$Position1,
              .data$Position2
            ),
            .data$Position2
          ),
          .data$Position2
        ),
        .data$Position2
      )
    ))

  problem_linked <- snp_df %>%
    dplyr::filter(.data$linked_snp_count == -1) %>%
    dplyr::mutate(Position = ifelse(
      .data$SNP1_effect < 0 & .data$SNP2_effect < 0,
      ifelse(
        .data$SNP1_effect != .data$SNP2_effect,
        ifelse(
          .data$SNP1_effect < .data$SNP2_effect,
          .data$Position1,
          .data$Position2
        ),
        .data$Position2
      ),
      ifelse(
        .data$SNP1_effect > 0 & .data$SNP2_effect > 0,
        ifelse(
          .data$SNP1_effect != .data$SNP2_effect,
          ifelse(
            .data$SNP1_effect > .data$SNP2_effect,
            .data$Position1,
            .data$Position2
          ),
          .data$Position2
        ),
        "problem"
      )
    ))

  snp_df <-
    rbind(no_linked, one_linked, problem_linked, more_linked) %>%
    dplyr::filter(Position != "problem") %>%
    dplyr::mutate(Position = as.integer(.data$Position)) %>%
    dplyr::mutate(window_start = .data$Position - 1000,
                  window_end = .data$Position + 1000) %>%
    dplyr::mutate(chr = as.character(.data$Locus1)) %>%
    dplyr::mutate(
      Effect = ifelse(
        .data$Position == .data$Position1,
        .data$SNP1_effect,
        .data$SNP2_effect
      ),
      P_value = ifelse(
        .data$Position == .data$Position1,
        .data$SNP1_pval,
        .data$SNP2_pval
      )
    ) %>%
    dplyr::select(
      .data$chr,
      .data$Position,
      .data$window_start,
      .data$window_end,
      .data$Effect,
      .data$P_value,
      .data$linked_snp_count
    )
  inner_join(snp_df, gff, by = "chr") %>%
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
    dplyr::mutate(
      Position = as.integer(.data$Position),
      Effect = as.numeric(.data$Effect),
      Locus1 = .data$chr
    ) %>%
    dplyr::select(-c(
      .data$chr,
      .data$start,
      .data$end,
      .data$window_start,
      .data$window_end
    ))
}

#' Parse SNP data and return genes
#'
#' @param all_data Merged association and effects data from merge_data()
#' @param LD Linkage disequilibrium data from parse_LD()
#' @param gff_file The path to a GFF file (accepts plain-text or gzipped plain-text)
#' @param window The search window for genes around the SNP
#' @param r_squared_cutoff The R^2 value used to determine SNP significance
#' @param num_cores The number of cores to use in parallelizing PAST
#' @importFrom rlang .data
#' @import dplyr
#' @import foreach
#' @import doParallel
#' @return A dataframe of genes from the SNP data
#' @export
#'
#' @examples
#' demo_association_file = system.file("extdata",
#' "association.txt.xz", package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata",
#' "effects.txt.xz", package = "PAST", mustWork = TRUE)
#' merged_data <- merge_data(demo_association_file, demo_effects_file)
#' demo_linkage_disequilibrium_file = system.file("extdata",
#' "LD.txt.xz", package = "PAST", mustWork = TRUE)
#' LD <- parse_LD(demo_linkage_disequilibrium_file)
#' demo_genes_file = system.file("extdata", "genes.gff.xz", package = "PAST", mustWork = TRUE)
#' genes <-parse_SNP(merged_data, LD, demo_genes_file, 1000, 0.8, 2)
parse_SNP <-
  function(all_data,
           LD,
           gff_file,
           window,
           r_squared_cutoff,
           num_cores) {
    full_gff <- read_gff(gff_file)
    chromosomes <- full_gff %>%
      dplyr::select(.data$chr) %>%
      dplyr::arrange(.data$chr) %>%
      unique()

    cl <- parallel::makeCluster(num_cores, outfile = "")
    registerDoParallel(cl)

    all_genes <- NULL

    # UP/DOWNSTREAM LOOP
    for (i in seq_along(c(1, 2))) {
      # print(ifelse(i == 1, "upstream", "dowstream"))

      # BEGIN PROCESSING BY CHROMOSOMES LOOP
      for (chromosome in names(LD)) {
        print(chromosome)

        if (chromosome %in% chromosomes$chr) {
          temp_data <-
            LD[[chromosome]] %>%
            dplyr::mutate(Marker1 = paste0("S", .data$Locus1, "_", .data$Position1)) %>%
            dplyr::mutate(Marker2 = paste0("S", .data$Locus1, "_", .data$Position2)) %>%
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

          # print("Parsing unlinked...")
          temp_data <-
            foreach(
              data = temp_data_list,
              .combine = rbind,
              .packages = c("dplyr", "PAST")
            ) %dopar% {
              parse_temp(data, r_squared_cutoff)
            }

          chr_data <-
            dplyr::filter(all_data, .data$Chr == as.integer(chromosome))

          # print("Getting effects...")
          # look up p-value and effect data for SNP1
          temp_data <-
            merge(temp_data, chr_data, by.x = "Marker1", by.y = "Marker") %>%
            dplyr::mutate(SNP1_pval = .data$p,
                          SNP1_effect = .data$Effect.x) %>%
            dplyr::select(
              .data$Locus1,
              .data$Position1,
              .data$Position2,
              .data$Site1,
              .data$Site2,
              .data$Dist_bp,
              .data$R.2,
              .data$Marker1,
              .data$SNP1_pval,
              .data$SNP1_effect,
              .data$Marker2
            )

          # look up p-value and effect data for SNP2
          temp_data <-
            merge(temp_data, chr_data, by.x = "Marker2", by.y = "Marker") %>%
            dplyr::mutate(SNP2_pval = .data$p,
                          SNP2_effect = .data$Effect.x) %>%
            dplyr::select(
              .data$Locus1,
              .data$Position1,
              .data$Position2,
              .data$Site1,
              .data$Site2,
              .data$Dist_bp,
              .data$R.2,
              .data$Marker1,
              .data$SNP1_pval,
              .data$SNP1_effect,
              .data$Marker2,
              .data$SNP2_pval,
              .data$SNP2_effect
            ) %>%
            dplyr::arrange(.data$Position1)

          index <- c(0, cumsum(abs(diff(
            temp_data$Site2
          )) > 1))
          temp_data_list <-
            split(temp_data, paste(temp_data$Position1, index))

          # print("Parsing chunks...")
          temp_data <-
            foreach(
              data = temp_data_list,
              .combine = rbind,
              .packages = c("dplyr", "PAST")
            ) %dopar% {
              parse_chunk(data, r_squared_cutoff)
            }

          split <- 1000
          snp_list <-
            split(temp_data,
                  rep(
                    seq(split),
                    length.out = nrow(temp_data),
                    each = ceiling(nrow(temp_data) / split)
                  ))

          # subset gff to only handle this chromosome
          gff <- dplyr::filter(full_gff, .data$chr == chromosome)

          # get genes in parallel
          # print("Finding genes...")
          chr_genes <-
            foreach(
              snp_chunk = snp_list,
              .combine = rbind,
              .packages = c("dplyr", "PAST")
            ) %dopar% {
              find_genes(gff, snp_chunk, window)
            }

          all_genes <- rbind(all_genes, chr_genes)
        }
      }
    }

    group_genes <- split(all_genes, f = all_genes$name)
    tagged_genes <-
      foreach(
        block = group_genes,
        .combine = rbind,
        .packages = c("dplyr", "PAST")
      ) %dopar% {
        select_gene_from_block(block)
      }
    
    stopCluster(cl)
    
    tagged_genes %>%
      dplyr::mutate(Chromosome = .data$Locus1, Gene = .data$name) %>%
      dplyr::select(.data$Chromosome,
                    .data$Position,
                    .data$Gene,
                    .data$Effect,
                    .data$P_value)
  }