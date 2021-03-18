#' Load GWAS data
#'
#' @param files a vector of files - either a single GWAS file, two files with
#' separate association and effects data with one marker per line, or two files
#' from the output of TASSEL
#' @param input_type one of "single", "two", or "TASSEL"
#' @param trait the name of the column containing the trait
#' @param marker the name of the column containing the marker
#' @param locus the name of the column containing the locus or chromosome
#' @param site the name of the column containing the site or position
#' @param p the name of the column containing the p-value
#' @param marker_R2 the name of the column containing the p-value
#' @param effect the name of the column containing the effect
#' @param effects_trait the name of the column containing the trait in the 
#' effects file; only used when input_type is "two" or "TASSEL"
#' @param effects_marker the name of the column containing the marker in the 
#' effects file; only used when input_type is "two" or "TASSEL"
#' @param effects_locus the name of the column containing the locus or 
#' chromosome in the effects file; only used when input_type is "two" or 
#' "TASSEL"
#' @param effects_site the name of the column containing the site or position
#' in the effects file; only used when input_type is "two" or "TASSEL"
#' @param update_progress used by the PAST Shiny app to report progress
#' @return GWAS data with a marker created from locus/chromosome and 
#' site/position, the original marker, locus/chromosome, site/position, p-value,
#' marker R^2, and effect for every marker in the input files
#' @export
#' @import data.table
#' @examples
#' association = system.file("extdata", "association.txt.gz", package = "PAST", 
#'   mustWork = TRUE)
#' effects = system.file("extdata", "effects.txt.gz", package = "PAST", 
#'   mustWork = TRUE)
#' single = system.file("extdata", "single_gwas.txt.gz", package = "PAST", 
#'   mustWork = TRUE)
#' effects_single = system.file("extdata", "effects.single_line.txt.gz", 
#'   package = "PAST", mustWork = TRUE)
#'   
#' one_file = c(single_file)
#' two_file = c(association_file, effects_single_file)
#' tassel = c(association_file, effects_file)
#' gwas_data <- load_GWAS_data(one_file, "single")
#' gwas_data <- load_GWAS_data(two_file, "two")
#' gwas_data <- load_GWAS_data(tassel, "TASSEL")
load_GWAS_data <- function(files,
                           input_type,
                           trait = "Trait",
                           marker = "Marker",
                           locus = "Locus",
                           site = "Site",
                           p = "p",
                           marker_R2 = "marker_R2",
                           effect = "Effect",
                           effects_trait = "Trait",
                           effects_marker = "Marker",
                           effects_locus = "Locus",
                           effects_site = "Site",
                           update_progress = NULL) {
  
  # Make a vector of column names provided by the user.
  association_columns <- c(trait, 
                          marker, 
                          locus, 
                          site, 
                          p, 
                          marker_R2)
  
  effects_columns <- c(effects_trait, 
                      effects_marker, 
                      effects_locus, 
                      effects_site, 
                      effect)
  
  # Compare the length of the input files vector and the type of analysis
  #   requested.
  # If length of the input files vector and input type match, continue.
  # Otherwise, fail with an error message.
  if (length(files) == 1 & input_type == "single") {
    
    # Read the file and select its columns.
    # Set the names to values used throughout PAST instead of what the user
    #   provided.
    # Create the marker column and drop the trait column.
    # Key the data by the chromosome and position columns and return.
    gwas <- data.table::fread(files[1], select = c(association_columns, effect))
    data.table::setnames(gwas, c(association_columns, effect), c("trait",
                                                                 "marker_original",
                                                                 "chromosome",
                                                                 "position",
                                                                 "p-value",
                                                                 "marker_R2",
                                                                 "effect"))
    gwas[,marker:=paste(chromosome, position, sep = "_")][,trait := NULL]
    data.table::setkeyv(gwas, c("chromosome", "position"))
    return(gwas)
    
  } else if (length(files) == 2 & input_type == "two") {
    
    # Read the file and select its columns.
    # Set the names to values used throughout PAST instead of what the user
    #   provided.
    # Create the marker column.
    # Key the data by the marker columns.
    associations <- data.table::fread(files[1], select = association_columns)
    data.table::setnames(associations, association_columns, c("trait",
                                                  "marker_original",
                                                  "chromosome",
                                                  "position",
                                                  "p-value",
                                                  "marker_R2"))
    associations[,marker:=paste(chromosome, position, sep = "_")]
    data.table::setkey(associations, marker)
    
    # Read the file and select its columns.
    # Set the names to values used throughout PAST instead of what the user
    #   provided.
    # Create the marker column and drop the trait, marker_original,
    #   chromosome, and position data.
    # Key the data by the marker columns.
    effects <- data.table::fread(files[2], select = effects_columns)
    data.table::setnames(effects, effects_columns, c("trait",
                                         "marker_original",
                                         "chromosome",
                                         "position",
                                         "effect"))
    effects[
      ,marker:=paste(chromosome, position, sep = "_")][
        ,c("trait", "marker_original", "chromosome", "position") := NULL
      ]
    data.table::setkey(effects, marker)
    
    # Merge associations with effects by the marker column and drop the 
    #   trait column.
    gwas <- associations[effects, on = .(marker)][,trait := NULL]
    
  } else if (length(files) == 2 & input_type == "TASSEL") {
    
    # Read the file and select its columns.
    # Set the names to values used throughout PAST instead of what the user
    #   provided.
    # Create the marker column.
    # Key the data by the marker columns.
    associations <- data.table::fread(files[1], select = association_columns)
    data.table::setnames(associations, association_columns, c("trait",
                                                  "marker_original",
                                                  "chromosome",
                                                  "position",
                                                  "p-value",
                                                  "marker_R2"))
    associations[,marker:=paste(chromosome, position, sep = "_")]
    data.table::setkey(associations, marker)
    
    # Read the file and select its columns.
    # Set the names to values used throughout PAST instead of what the user
    #   provided.
    # Create the marker column.
    # Find out which markers are biallelic.
    # Reshape the data so that both effects for each marker on on the same
    #   line, only keeping the marker column and the first effect.
    #   (TASSEL's secondary effects are 0.)
    # Rename the automatically named "1" column to "effect".
    # Key the data by the marker columns.
    effects <- data.table::fread(files[2], select = effects_columns)
    data.table::setnames(effects, effects_columns, c("trait",
                                         "marker_original",
                                         "chromosome",
                                         "position",
                                         "effect"))
    effects[,marker:=paste(chromosome, position, sep = "_")]
    biallelic <- effects[, .(.N), by = .(marker)][N == 2]
    effects <- data.table::dcast(effects, 
                                marker ~ rowid(marker), 
                                value.var = c("effect"))[
                                  biallelic, c("marker", "1")
                                  ]
    data.table::setnames(effects, "1", "effect")
    data.table::setkey(effects, marker)
    
    # Merge associations with effects by the marker column and drop the 
    #   trait column.
    gwas <- associations[effects, on = .(marker)][,trait := NULL]
    
  } else {
    stop("Length of files and input type are incorrect.")
  }
  # Key the data by the chromosome and position columns and return.
  data.table::setkeyv(gwas, c("chromosome", "position"))
  return(gwas)
}
