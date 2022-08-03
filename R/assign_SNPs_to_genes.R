#' Assign SNPs to genes
#'
#' @param gwas_data GWAS data from load_GWAS_data()
#' @param LD linkage disequilibrium data from load_LD()
#' @param GFF_file the path to a GFF-formatted of annotations.
#' @param filter_type keep only this type of feature from the GFF
#' @param window the search window for genes around the SNP
#' @param r_squared_cutoff the R^2 value used to determine SNP significance
#' @param name_attribute the field in the 9th column of the GFF file that
#' contains the desired gene name (for example, the 9th column might look like
#' ID=NNN;Name=NNN;description=NNN - to use the name in the Name field, provide
#' "Name")
#' @param debug return all SNPs associated with a gene
#' @import data.table
#' @import stringr
#' @return a data.table of SNP-gene assignments, reduced to the SNP that best
#' represents all SNPs assigned to the same gene for each gene
#' @export
#'
#' @examples
#' example("load_GWAS_data")
#' example("load_LD")
#' GFF = system.file("extdata", "genes.gff.gz", package = "PAST2",
#'   mustWork = TRUE)
#' genes <- assign_SNPs_to_genes(gwas_data, LD, GFF, "gene", 1000, 0.8, "Name")
assign_SNPs_to_genes <- function(
    gwas_data,
    LD,
    GFF_file,
    filter_type = "gene",
    window,
    r_squared_cutoff,
    name_attribute,
    debug = FALSE) {
  
  # Assign some variables to NULL to pass checks.
  type <- name <- seqid <- locus <- marker1 <- position1 <- marker2 <- NULL
  position2 <- site1 <- represents <- distance <- r_squared <- site2 <- NULL
  difference <- block <- linked_snp_count <- effect2 <- SNP <- effect1 <- NULL
  position <- `p-value` <- `p-value1` <- effect <- `p-value2` <- NULL
  chromosome <- marker <- start <- end <- .  <- temp_position <- NULL
  temp_site <- NULL
  
  # PAST technically accepts a data.table of genes or a path to a GFF file.
  # However, PAST only accepts a data.table because doing so makes working
  #   with PAST Shiny easier. Users should pass a path.
  # PAST only needs five GFF columns.
  # Read the GFF-file, selecting only those columns and naming them.
  # Get the names of all fields in the attributes column.
  # Split the attributes column by ";", naming the columns based on the
  #   attribute name.
  # Set the name of the column containing the gene name to "name".
  # Drop the other attribute columns.
  # Add the search window to the gene coordinates.
  # Key the GFF data.table by sequence ID, start, and end in preparation for
  #   finding intersections with the tag SNPs
  if (data.table::is.data.table(GFF_file) == FALSE) {
    GFF_columns <- c("seqid", "type", "start", "end", "attributes")
    GFF <- data.table::as.data.table(
      read.table(
        GFF_file,
        comment.char = "#",
        sep = "\t",
        header = FALSE,
        quote=""
      )
    )[, c("V1", "V3", "V4", "V5", "V9")]
    data.table::setnames(GFF, c("V1", "V3", "V4", "V5", "V9"), GFF_columns)
    GFF <- GFF[type == filter_type]
    attributes = GFF[, data.table::tstrsplit(attributes, ";", fixed = TRUE)]
    names = apply(
      attributes,
      1,
      function(row) {
        names(attributes)[
          which(
            stringr::str_detect(
              row,
              paste0(name_attribute, "=")
            )
          )
        ]
      }
    )
    
    for (row in 1:nrow(GFF)){
      data.table::set(
        GFF,
        i = row,
        j = "name",
        value = stringr::str_replace(
          as.character(attributes[row, names[row], with = FALSE]),
          paste0(name_attribute, "="),
          "")
      )
    }
    GFF[, attributes := NULL]
    GFF[, seqid := as.character(seqid)]
    GFF[, start := start - window]
    GFF[, end := end + window]
    data.table::setkey(GFF, seqid, start, end)
  }
  
  genes <- data.table::rbindlist(
    lapply(seq_along(c(1, 2)), function(direction) {
      # Subset the LD data by discarding chromosomes not in the GFF.
      initial_data <- LD[locus %in% GFF$seqid]
      if (direction == 1) {
        initial_data[, marker1 := paste(locus, position1, sep = "_")]
        initial_data[, marker2 := paste(locus, position2, sep = "_")]
      } else {
        initial_data[, marker1 := paste(locus, position2, sep = "_")]
        initial_data[, marker2 := paste(locus, position1, sep = "_")]
        initial_data[
          ,
          c("temp_position", "temp_site") := .(position1, site1)
        ]
        initial_data[
          ,
          c("position1", "site1") := .(position2, site2)
        ]
        initial_data[
          , c("position2", "site2") := .(temp_position, temp_site)
        ]
        initial_data[
          ,
          c("temp_position", "temp_site") := NULL
        ]
      }
      data.table::setorder(initial_data, locus, position1, site1)
      
      # If debug is set, keep which SNPs are associated with the linkage.
      if (debug) {
        initial_data[, represents := list(character(0))]
        initial_data[, represents := {
          list(list(c(unique(.SD$marker1), unique(.SD$marker2))))
        }, by = .(locus, position1)]
      }
      
      # Group the data by locus and position1.
      # For every group, check to see whether all of the rows are
      #   unlinked.
      # If all rows are unlinked, find the row with the smallest distance
      #   between SNP1 and SNP2.
      # Otherwise, drop the unlinked rows.
      linkage <- initial_data[, {
        SNPs <- c(unique(.SD$marker1), unique(.SD$marker2))
        if (sum(.SD$r_squared <= r_squared_cutoff) == .N) {
          .SD[which.min(distance)]
        } else {
          .SD[r_squared > r_squared_cutoff]
        }
      }, by = .(locus, position1)]
      
      # Merge GWAS data markers, effects, and p-values with linkage data
      #   marker in the GWAS data is equal to marker1 in the linkage
      #   data.
      # Set the names to show that this data is associated with SNP1.
      # Drop markers where there was no match.
      linkage_with_effects <- gwas_data[
        ,
        c("marker", "effect", "p-value")
      ][linkage, on = .(marker = marker1)]
      data.table::setnames(
        linkage_with_effects,
        c("marker", "effect", "p-value"),
        c("marker1", "effect1", "p-value1")
      )
      linkage_with_effects <- stats::na.omit(
        linkage_with_effects,
        cols = "effect1")
      
      # Merge GWAS data markers, effects, and p-values with linkage data
      #   where marker in the GWAS data is equal to marker2 in the linkage
      #   data.
      # Set the names to show that this data is associated with SNP2.
      # Drop markers where there was no match.
      linkage_with_effects <- gwas_data[
        ,
        c("marker", "effect", "p-value")
      ][linkage_with_effects, on = .(marker = marker2)]
      data.table::setnames(
        linkage_with_effects,
        c("marker", "effect", "p-value"),
        c("marker2", "effect2", "p-value2")
      )
      data.table::setcolorder(
        linkage_with_effects,
        c(
          "locus",
          "marker1",
          "position1",
          "site1",
          "effect1",
          "p-value1",
          "marker2",
          "position2",
          "site2",
          "effect2",
          "p-value2",
          "distance",
          "r_squared"
        )
      )
      data.table::setorder(linkage_with_effects, locus, position1)
      
      # Mark blocks of SNPs where position1 is the same and site2 is
      #   consecutive.
      # Select only blocks of SNPs where there is more than one SNP with
      #   the same marker1.
      # Order these blocks by locus, position1, site1, and site2.
      # Omit any SNPs that have no effect2 data.
      # Create a difference column to detect consecutive site2s.
      # Set the first row's difference to 0.
      # Create a block column based on the cumulative sum of the
      #   difference column.
      # Set the number of linked SNPs. SNPs with R^2 less than the cutoff
      #   are unlinked. SNPs with R^2 than the cutoff are considered
      #   linked to the number of SNPs in the block.
      # Blocks with an equal number of positive and negative SNPs are
      #   marked with a number of linked SNPs equal to -1.
      block_SNPs <- linkage_with_effects[, SNPs := .N, by = .(marker1)][
        SNPs > 1
      ]
      data.table::setorder(block_SNPs, locus, position1, site1, site2)
      block_SNPs <- stats::na.omit(block_SNPs, cols = "effect2")
      block_SNPs[
        ,
        difference := abs(
          site2 - data.table::shift(
            site2,
            1,
            "lag",
            fill = 0
          )
        ) > 1 |
          position1 != data.table::shift(
            position1,
            1,
            "lag",
            fill = 0
          )
      ]
      block_SNPs[1L, difference := 0]
      block_SNPs[, block := cumsum(difference)]
      block_SNPs[, linked_snp_count := -999]
      
      # Process the blocks to find the single SNP that best represents
      #   each block.
      # Group by the locus, position1, and block columns.
      # Count negative and positive effects in the block.
      # If there are more positive effects, sort the block by descending
      #   effect2 and distance.
      # If there are more negative effects, sort the block by effect2 and
      #   distance.
      # If there are an equal number of positive and negative effects,
      #   then check the p-value of the first SNP in the block.
      # If that p-value is positive, sort the block by descending effect2
      #   and distance.
      # If that p-value is negative, sort the block by effect2 and
      #   distance.
      # Return the first row of the block, which is the SNP with the
      #   largest most negative effect if there were more negative SNPs or
      #   the most positive effect if there were more positive SNPs.
      block_SNPs <- block_SNPs[, {
        SNPs_in_block <- .N
        negative <- .SD[effect2 < 0][, .N]
        positive <- .SD[effect2 > 0][, .N]
        subdata <- data.table::copy(.SD)
        subdata[r_squared < r_squared_cutoff,
                linked_snp_count := as.numeric(0)]
        subdata[r_squared >= r_squared_cutoff,
                linked_snp_count := as.numeric(SNPs_in_block)]
        if (positive > negative) {
          data.table::setorder(subdata, -effect2, distance)
        } else if (negative > positive) {
          data.table::setorder(subdata, effect2, distance)
        } else if (negative == positive) {
          subdata[, linked_snp_count := as.numeric(-1)]
          if (subdata[1, effect1] > 0) {
            data.table::setorder(subdata, -effect2, distance)
          } else {
            data.table::setorder(subdata, effect2, distance)
          }
        }
        subdata[1L]
      },
      by = block]
      block_SNPs[, c("block", "SNPs", "difference") := NULL]
      data.table::setcolorder(
        block_SNPs,
        c(
          "locus",
          "marker1",
          "position1",
          "site1",
          "effect1",
          "p-value1",
          "marker2",
          "position2",
          "site2",
          "effect2",
          "p-value2",
          "distance",
          "r_squared",
          "linked_snp_count"
        )
      )
      data.table::setorder(block_SNPs, locus, position1, site1, site2)
      
      # Determine whether single SNPs are linked to their SNP2 or not.
      # Select only blocks of SNPs where there is only one SNP at that
      #   position.
      # Set the number of linked SNPs. SNPs with R^2 less than the cutoff
      #   are unlinked. SNPs with R^2 than the cutoff are considered
      #   linked.
      single_SNPs <- linkage_with_effects[, SNPs := .N, by = .(marker1)][
        SNPs == 1
      ][r_squared < r_squared_cutoff, linked_snp_count := 0][
        r_squared >= r_squared_cutoff, linked_snp_count := 1]
      single_SNPs[, SNPs := NULL]
      data.table::setcolorder(
        single_SNPs,
        c(
          "locus",
          "marker1",
          "position1",
          "site1",
          "effect1",
          "p-value1",
          "marker2",
          "position2",
          "site2",
          "effect2",
          "p-value2",
          "distance",
          "r_squared",
          "linked_snp_count"
        )
      )
      
      # Bind the SNPs representing blocks of SNPs to the single and
      #   unlinked SNPs.
      # Order the data.
      representative_SNPs <- data.table::rbindlist(
        list(
          block_SNPs,
          single_SNPs
        )
      )
      data.table::setorder(
        representative_SNPs,
        locus,
        position1,
        site1,
        site2
      )
      
      # Find tag SNPs - either SNP1 or SNP2 which represents the pair.
      # If effect1 and effect2 are equal, always select SNP2.
      tagSNPs <- data.table::copy(representative_SNPs)[effect1 == effect2, SNP := "SNP2"]
      
      # If SNP1 and SNP2 are not linked, select SNP1.
      tagSNPs[linked_snp_count == 0 & is.na(SNP), SNP := "SNP1"]
      
      # If SNP1 and SNP2 are linked to each other and their effects have
      #   the same sign, choose the SNP with the highest absolute effect.
      tagSNPs[linked_snp_count == 1 &
                effect1 * effect2 > 0 &
                abs(effect1) > abs(effect2) &
                is.na(SNP),
              SNP := "SNP1"]
      tagSNPs[linked_snp_count == 1 &
                effect1 * effect2 > 0 &
                abs(effect1) < abs(effect2) &
                is.na(SNP),
              SNP := "SNP2"]
      
      # If SNP1 and SNP2 are linked to each other and their effects have
      #   different signs, choose the SNP with the lower p-value.
      tagSNPs[linked_snp_count == 1 &
                effect1 * effect2 < 0 &
                `p-value1` > `p-value2` &
                is.na(SNP),
              SNP := "SNP2"]
      tagSNPs[linked_snp_count == 1 &
                effect1 * effect2 < 0 &
                `p-value1` < `p-value2` &
                is.na(SNP),
              SNP := "SNP1"]
      
      # If SNP1 and SNP2 are linked to each other and their effects have
      #   different signs and they have equal p-values, choose the SNP
      #   with the largest absolute effect.
      tagSNPs[linked_snp_count == 1 &
                effect1 * effect2 < 0 &
                `p-value1` == `p-value2` &
                abs(effect1) > abs(effect2) &
                is.na(SNP),
              SNP := "SNP1"]
      tagSNPs[linked_snp_count == 1 &
                effect1 * effect2 < 0 &
                `p-value1` == `p-value2` &
                abs(effect1) < abs(effect2) &
                is.na(SNP),
              SNP := "SNP2"]
      
      # If the linkage represents a block of SNPs and and their effects
      #   have different  signs or the distance between the SNPs is
      #   greater than the search window, choose SNP2.
      tagSNPs[linked_snp_count > 1 & effect1 * effect2 < 0 & is.na(SNP), SNP := "SNP2"]
      tagSNPs[linked_snp_count > 1 & distance >= window & is.na(SNP), SNP := "SNP2"]
      
      # If the linkage represents a block of SNPs and their effects have
      #   the same signs and the distance between the SNPs is smaller than
      #   the search window, choose the SNP with the highest absolute
      #   effect.
      tagSNPs[linked_snp_count > 1 &
                effect1 * effect2 > 0 &
                distance < window &
                abs(effect1) > abs(effect2) &
                is.na(SNP),
              SNP := "SNP1"]
      tagSNPs[linked_snp_count > 1 &
                effect1 * effect2 > 0 &
                distance < window &
                abs(effect1) < abs(effect2) &
                is.na(SNP),
              SNP := "SNP2"]
      
      # If the linkage represents a block of SNPs with the same number of
      #   positive and negative effects and their effects have the same
      #   sign, choose the SNP with the highest absolute effect.
      tagSNPs[linked_snp_count == -1 &
                effect1 * effect2 > 0 &
                abs(effect1) > abs(effect2) &
                is.na(SNP),
              SNP := "SNP1"]
      tagSNPs[linked_snp_count == -1 &
                effect1 * effect2 > 0 &
                abs(effect1) < abs(effect2) &
                is.na(SNP),
              SNP := "SNP2"]
      
      # If the linkage represents a block of SNPs with the same number of
      #   positive and negative effects and their effects have different
      #   signs, note that the linkage is problematic.
      tagSNPs[linked_snp_count == -1 &
                effect1 * effect2 < 0 &
                is.na(SNP),
              SNP := "PROBLEM"]
      
      # Set position, effect, and p-value to the values for the tag SNP.
      tagSNPs[SNP == "SNP1", position := position1]
      tagSNPs[SNP == "SNP1", effect := effect1]
      tagSNPs[SNP == "SNP1", `p-value` := `p-value1`]
      tagSNPs[SNP == "SNP2", position := position2]
      tagSNPs[SNP == "SNP2", effect := effect2]
      tagSNPs[SNP == "SNP2", `p-value` := `p-value2`]
      
      # Change the name of the locus column to chromosome, drop rows with
      #   no position, and select only needed columns.
      # Add start and end columns equal to the position column in
      #   preparation for finding intersections with the GFF.
      # Order the data.
      data.table::setnames(tagSNPs, "locus", "chromosome")
      tagSNPs <- stats::na.omit(tagSNPs, cols = "position")
      if (debug) {
        tagSNPs <- tagSNPs[
          ,
          .(
            chromosome,
            position,
            effect,
            `p-value`,
            distance,
            linked_snp_count,
            represents
          )
        ]
      } else {
        tagSNPs <- tagSNPs[
          ,
          .(
            chromosome,
            position,
            effect,
            `p-value`,
            distance,
            linked_snp_count
          )
        ]
      }
      tagSNPs[, start := position]
      tagSNPs[, end := position]
      data.table::setorder(
        tagSNPs,
        chromosome,
        position,
        distance,
        linked_snp_count
      )
      
      # Key the tag SNP data.table by sequence ID, start, and end in
      #   preparation for finding intersections with the GFF.
      data.table::setkey(tagSNPs, chromosome, start, end)
      
      # Find overlaps between the tag SNPs and GFF, dropping rows from the
      #   tag SNPs with no match in the GFF.
      # Only keep specified columns.
      # Order the data, and set number of linked SNPs equal to 1 for any
      #   SNPs with a number of linked SNPs less than or equal to 0.
      if (debug) {
        genes <- data.table::foverlaps(
          tagSNPs,
          GFF,
          nomatch = NULL)[
            ,
            .(
              chromosome,
              position,
              effect,
              `p-value`,
              linked_snp_count,
              name,
              represents
            )
          ]
      } else {
        genes <- data.table::foverlaps(
          tagSNPs,
          GFF,
          nomatch = NULL)[
            ,
            .(
              chromosome,
              position,
              effect,
              `p-value`,
              linked_snp_count,
              name
            )
          ]
      }
      
      data.table::setorder(
        genes,
        chromosome,
        position,
        effect,
        `p-value`,
        linked_snp_count,
        name
      )
      genes[linked_snp_count <= 0, linked_snp_count := 1]
      genes
    }
    ))
  
  genes
  # Find the SNP-gene pairing that best represents all SNPs linked to the same
  #   gene.
  # Sum the number of linked SNPs with negative effects and the number of
  #   linked SNPs with positive effects.
  # If there are more positive SNPs linked to the gene, sort by descending
  #   effect and ascending p-value. Select the first SNP in the
  #   sub-data.table.
  # If there are more negative SNPs linked to the gene, sort by ascending
  #   effect and p-value. Select the first SNP in the sub-data.table.
  # If there are an equal number of SNPs with positive and negative effects
  #   linked to the gene, select the SNP with the largest absolute effect.
  # Set the number of linked SNPs to all SNPs linked to the gene.
  # Set the column order and order the data.
  genes <- genes[, {
    negative <- .SD[effect < 0, sum(linked_snp_count)]
    positive <- .SD[effect > 0, sum(linked_snp_count)]
    subdata <- data.table::copy(.SD)
    if (positive > negative) {
      data.table::setorder(subdata, -effect, `p-value`, -linked_snp_count, position)
      subdata <- subdata[1L]
    } else if (negative > positive) {
      subdata = data.table::setorder(subdata, effect, `p-value`, -linked_snp_count, position)
      sub2 = subdata %>% arrange(effect, `p-value`, -linked_snp_count, position)
      subdata <- subdata[1L]
    } else if (negative == positive) {
      max_positive <- data.table::setorder(
        subdata[effect >= 0],
        -effect,
        `p-value`,
        -linked_snp_count, 
        position
      )[1L]
      max_negative <- data.table::setorder(
        subdata[effect <= 0],
        effect,
        `p-value`,
        -linked_snp_count,
        position
      )[1L]
      
      if (max_positive[, effect] >= abs(max_negative[, effect])) {
        subdata <- max_positive
      } else {
        subdata <- max_negative
      }
    }
    if (debug) {
      SNPs <- unique(unlist(.SD$represents))
      subdata[, c("linked_snp_count", "represents") :=
                list(negative + positive, SNPs)]
    } else {
      subdata[, linked_snp_count := negative + positive]
    }
    
  },
  by = name]
  
  data.table::setcolorder(
    genes,
    c(
      "chromosome",
      "position",
      "effect",
      "p-value",
      "linked_snp_count",
      "name"
    )
  )
  data.table::setorder(
    genes,
    chromosome,
    position,
    effect,
    `p-value`,
    linked_snp_count,
    name
  )
  
  # Retrieve the original marker and rename the marker_original column to
  #   marker.
  genes <- gwas_data[, c("marker", "marker_original")][
    genes[, marker := paste(chromosome, position, sep = "_")],
    on = .(marker = marker)
  ]
  genes[, marker := NULL]
  data.table::setnames(genes, "marker_original", "marker")
  genes
}
