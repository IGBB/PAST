# These aren't really global, but R CMD check
# doesn't understand how to deal with variable
# created by foreach loops
globalVariables("snp_chunk")
globalVariables("chunk")

#' Determine Linkage
#'
#' @param chunk A chunk of data to be processed
#' @param r_squared_cutoff The R^2 value to check against
#' @importFrom rlang .data
#' @import dplyr
#' @return Either the first unlinked SNP or a set of linked SNPs
determine_linkage <- function(chunk, r_squared_cutoff) {
    lapply(chunk, function(data) {
        data <- data %>% dplyr::arrange(.data$Dist_bp)
        linked <- data %>% dplyr::filter(.data$R.2 >= r_squared_cutoff)
        unlinked <- data %>% dplyr::filter(.data$R.2 < r_squared_cutoff)
        if (nrow(unlinked) == nrow(data)) {
            unlinked[1, ]
        } else {
            linked %>% filter(.data$R.2 > r_squared_cutoff)
        }
    }) %>% dplyr::bind_rows()
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
                dplyr::arrange(desc(.data$SNP2_effect), .data$Dist_bp)
            
            # get top row
            chunk <- chunk[1, ] %>%
                dplyr::mutate(linked_snp_count = -1)
        } else{
            # sort in ascending order
            chunk <- chunk %>%
                dplyr::arrange(.data$SNP2_effect, .data$Dist_bp)
            
            # get top row
            chunk <- chunk[1, ] %>%
                dplyr::mutate(linked_snp_count = -1)
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
#' @import GenomicRanges
#' @import S4Vectors
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
                                                     .data$largest_effect),
                      pvalue_equal = ifelse(.data$SNP1_pvalue == .data$SNP2_pvalue,
                                            TRUE,
                                            FALSE),
                      pvalue_not_equal_use = ifelse(.data$SNP1_pvalue > .data$SNP2_pvalue,
                                                    "SNP2",
                                                    "SNP1"),
                      single_pvalue_use = ifelse(.data$pvalue_equal == FALSE,
                                                 .data$pvalue_not_equal_use,
                                                 "PROBLEM"),
                      single_use = ifelse(.data$effects_same_sign == TRUE,
                                          .data$effects_same_sign_use,
                                          .data$pvalue_not_equal_use),
                      distance_use = ifelse(.data$Dist_bp < window,
                                            .data$effects_same_sign_use,
                                            "SNP2"),
                      many_use = ifelse(.data$effects_same_sign == TRUE,
                                        .data$distance_use,
                                        "SNP2"),
                      problem_use = ifelse(.data$effects_same_sign == TRUE,
                                           .data$effects_same_sign_use,
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
                      .data$linked_snp_count,
                      .data$Marker_original
        )
    
    # handle unlinked SNPs
    conditions <- rbind(conditions %>% dplyr::filter(.data$unlinked == FALSE),
                        conditions %>% dplyr::filter(.data$unlinked == TRUE) %>%
                            dplyr::mutate(position = .data$Position1,
                                          effect = .data$SNP1_effect,
                                          p.value = .data$SNP1_pvalue)
    ) %>%
        dplyr::select(-.data$unlinked)
    
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
        dplyr::select(-.data$single_linked,
                      -.data$single_use)
    
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
        dplyr::select(-.data$many_linked,
                      -.data$many_use)
    
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
        dplyr::mutate(marker = .data$Marker_original) %>%
        dplyr::select(.data$chromosome,
                      .data$position,
                      .data$marker,
                      .data$effect,
                      .data$p.value,
                      .data$distance,
                      .data$linked_snp_count) %>%
        dplyr::mutate(start = .data$position - window,
                      stop = .data$position + window) %>%
        dplyr::arrange(.data$position) %>% filter(!(is.na(.data$position)))
    
    # assign genes
    gr1 = with(tagSNPs, GRanges(chromosome, IRanges(start = position, 
                                                    end = position)))
    gr2 = with(gff, GRanges(seqid, IRanges(start = start-window, 
                                           end = end+window, 
                                           names = Name)))
    overlaps = GenomicRanges::findOverlaps(query = gr1, subject = gr2)
    data.frame(tagSNPs[S4Vectors::queryHits(overlaps),], 
               name=gff[S4Vectors::subjectHits(overlaps),"Name"]) %>%
        dplyr::select(-c(
            .data$distance,
            .data$start,
            .data$stop
        ))
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
#' example("load_GWAS_data")
#' example("load_LD")
#' demo_genes_file = system.file("extdata", "genes.gff",
#'   package = "PAST", mustWork = TRUE)
#' filter_type = c("gene")
#' genes <-assign_SNPs_to_genes(gwas_data, LD, demo_genes_file, filter_type, 1000, 0.8, 2)
assign_SNPs_to_genes <-
    function(gwas_data,
             LD,
             gff_file,
             filter_type,
             window,
             r_squared_cutoff,
             num_cores) {
        
        full_gff <- as.data.frame(rtracklayer::readGFF(gff_file, filter=list(type=filter_type))) %>%
            select(.data$seqid,
                   .data$start,
                   .data$end,
                   .data$Name)
        chromosomes <- as.character.factor(
            full_gff %>%
                dplyr::select(.data$seqid) %>%
                dplyr::arrange(.data$seqid) %>%
                unique() %>%
                dplyr::pull(.data$seqid)
        )
        
        cl <- parallel::makeCluster(num_cores)
        clusterEvalQ(cl, {library(dplyr); library(GenomicRanges)})

        all_genes <- NULL
        
        # UP/DOWNSTREAM LOOP
        for (i in seq_along(c(1, 2))) {
            
            # BEGIN PROCESSING BY CHROMOSOMES LOOP
            for (chromosome in names(LD)) {
                
                if (chromosome %in% chromosomes) {
                    temp_data <-
                        LD[[chromosome]] %>%
                        dplyr::mutate(Marker1 = paste0(.data$Locus,
                                                       "_",
                                                       .data$Position1),
                                      Marker2 = paste0(.data$Locus,
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
                    split <- 1000
                    temp_data_list <- split(temp_data, temp_data$Position1)
                    temp_data_list_split <- split(temp_data_list,
                                                  ceiling(seq_along(temp_data_list)
                                                          /split))
                    
                    temp_data <- parLapply(cl, 
                                           temp_data_list_split, 
                                           determine_linkage, 
                                           r_squared_cutoff = r_squared_cutoff) %>%
                        bind_rows()
                    
                    gwas_data_for_chromosome <-
                        dplyr::filter(gwas_data,
                                      .data$Chr == as.integer(chromosome))
                    
                    # look up p-value and effect data for SNP1
                    temp_data <-
                        merge(gwas_data_for_chromosome, temp_data,
                              by.x = "Marker",
                              by.y = "Marker1") %>%
                        dplyr::mutate(Marker1 = .data$Marker,
                                      SNP1_pvalue = .data$p,
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
                            .data$Marker2,
                            .data$Marker_original
                        )
                    
                    # look up p-value and effect data for SNP2
                    temp_data <-
                        merge(gwas_data_for_chromosome, temp_data,
                              by.x = "Marker",
                              by.y = "Marker2",
                              all.y = TRUE) %>%
                        dplyr::mutate(Marker2 = .data$Marker,
                                      SNP2_pvalue = .data$p,
                                      SNP2_effect = .data$Effect.x,
                                      Marker_original = .data$Marker_original.x) %>%
                        dplyr::select(
                            .data$Locus,
                            .data$Position1,
                            .data$Site1,
                            .data$Position2,
                            .data$Site2,
                            .data$Dist_bp,
                            .data$R.2,
                            .data$Marker1,
                            .data$SNP1_pvalue,
                            .data$SNP1_effect,
                            .data$Marker2,
                            .data$SNP2_pvalue,
                            .data$SNP2_effect,
                            .data$Marker_original
                        ) %>%
                        dplyr::arrange(.data$Position1)
                    
                    singles = temp_data %>% dplyr::group_by(.data$Marker1) %>%
                        dplyr::summarise(count = n()) %>%
                        dplyr::filter(count == 1)
                    
                    linked_to_one <- temp_data %>% 
                        filter(.data$Marker1 %in% singles$Marker1, 
                               .data$R.2 > r_squared_cutoff) %>%
                        mutate(linked_snp_count = 1)
                    
                    linked_to_none <- temp_data %>% 
                        filter(.data$Marker1 %in% singles$Marker1, 
                               .data$R.2 <= r_squared_cutoff) %>%
                        mutate(linked_snp_count = 0)
                    
                    blocks <- temp_data %>% 
                        filter(!(.data$Marker1 %in% singles$Marker1))
                    
                    blocks <- blocks %>% filter(!(is.na(.data$SNP2_effect)))
                    
                    index <- c(0, cumsum(abs(diff(blocks$Site2)) > 1))
                    temp_data_list <- split(blocks, 
                                            paste(blocks$Position1, index))
                    
                    temp_data <- parLapply(cl, 
                                           temp_data_list, 
                                           find_representative_SNP, 
                                           r_squared_cutoff = r_squared_cutoff) %>%
                        bind_rows()
                    
                    temp_data <- rbind(temp_data, linked_to_one, linked_to_none) %>% 
                        arrange(.data$Position1)
                    
                    split <- 4
                    snp_list <-
                        split(temp_data,
                              rep(seq_len(split),
                                  length.out = nrow(temp_data),
                                  each = ceiling(nrow(temp_data) / split)
                              ))
                    
                    # subset gff to only handle this chromosome
                    gff <- dplyr::filter(full_gff, .data$seqid == chromosome) %>%
                        mutate(seqid = as.character.factor(.data$seqid))
                    
                    # get genes in parallel
                    chromosome_genes <- parLapply(cl, 
                                                  snp_list, 
                                                  assign_chunk, 
                                                  gff = gff, 
                                                  window = window) %>%
                        bind_rows()
                    
                    all_genes <- rbind(all_genes, 
                                       chromosome_genes %>% 
                                           arrange(.data$position))
                }
            }
        }
        
        group_by_gene <- split(all_genes, f = all_genes$name)
        representative_genes <-parLapply(cl, 
                                         group_by_gene, 
                                         find_representative_SNP_gene_pairing) %>%
            bind_rows()
        stopCluster(cl)
        representative_genes
    }
