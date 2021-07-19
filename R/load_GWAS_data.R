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
                           mode = "homozygous",
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
                           effects_allele = "Allele",
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
                       effects_allele,
                       effect)
  
  # Compare the length of the input files vector and the type of analysis
  #   requested.
  # If length of the input files vector and input type match, continue.
  # Otherwise, fail with an error message.
  if (length(files) == 1 & input_type == "single") {
    
    # Check the headers of the data to be sure that all user-requested
    #   columns exist.
    # If they don't, throw an error and exit.
    header = data.table::fread(files[1], nrows = 1, header = FALSE)
    column_check <- c(association_columns, effect) %in% unlist(header)
    names(column_check) = c(association_columns, effect)
    if (!all(column_check)) {
      stop(paste0("Could not find the following columns in data: ", 
                  paste(unlist(attr(column_check[column_check == FALSE], "names")), collapse = ", "),
                  "\n",
                  "Column names in data are: ",
                  paste(unlist(header), collapse = ", "))
      )
    }
    
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
    
    # Check the headers of the data to be sure that all user-requested
    #   columns exist.
    # If they don't, throw an error and exit.
    association_header = data.table::fread(files[1], nrows = 1, header = FALSE)
    column_check <- association_columns %in% unlist(association_header)
    names(column_check) = association_columns
    if (!all(column_check)) {
      stop(paste0("Could not find the following columns in association data: ", 
                  paste(unlist(attr(column_check[column_check == FALSE], "names")), collapse = ", "),
                  "\n",
                  "Column names in association data are: ",
                  paste(unlist(association_header), collapse = ", "))
      )
    }
    
    effects_header = data.table::fread(files[2], nrows = 1, header = FALSE)
    column_check <- effects_columns %in% unlist(effects_header)
    names(column_check) = effects_columns
    if (!all(column_check)) {
      stop(paste0("Could not find the following columns in effects data: ", 
                  paste(unlist(attr(column_check[column_check == FALSE], "names")), collapse = ", "),
                  "\n",
                  "Column names in effects data are: ",
                  paste(unlist(effects_header), collapse = ", "))
      )
    }
    
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
    
    # Check the headers of the data to be sure that all user-requested
    #   columns exist.
    # If they don't, throw an error and exit.
    association_header = data.table::fread(files[1], nrows = 1, header = FALSE)
    column_check <- association_columns %in% unlist(association_header)
    names(column_check) = association_columns
    if (!all(column_check)) {
      stop(paste0("Could not find the following columns in association data: ", 
                  paste(unlist(attr(column_check[column_check == FALSE], "names")), collapse = ", "),
                  "\n",
                  "Column names in association data are: ",
                  paste(unlist(association_header), collapse = ", "))
      )
    }
    
    effects_header = data.table::fread(files[2], nrows = 1, header = FALSE)
    column_check <- effects_columns %in% unlist(effects_header)
    names(column_check) = effects_columns
    if (!all(column_check)) {
      stop(paste0("Could not find the following columns in effects data: ", 
                  paste(unlist(attr(column_check[column_check == FALSE], "names")), collapse = ", "),
                  "\n",
                  "Column names in effects data are: ",
                  paste(unlist(effects_header), collapse = ", "))
      )
    }
    
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
    # Count the number of occurrences of each marker.
    # Keep all homozygous data and all heterozygous data.
    #   Heterozygous data has three occurrences per marker and 
    #   is noted by specific IUPAC codes.
    # Key the data by the marker columns.
    effects <- data.table::fread(files[2], select = effects_columns)
    data.table::setnames(effects, effects_columns, c("trait",
                                                     "marker_original",
                                                     "chromosome",
                                                     "position",
                                                     "allele",
                                                     "effect"))
    effects[,marker:=paste(chromosome, position, sep = "_")]
    effects[, num_occurrences := .(.N), by = .(marker)]
    effects = effects[order(effect), .SD, marker]

    # Inform user about mode.
    if (mode == "homozygous") {
      message("homozygous mode: retain only biallelic homozygous data")
    } else if (mode == "decreasing") {
      message("decreasing mode: retain marker with largest negative effect for heterozygous data")
      one = 1
      row_choice = "one"
    } else if (mode == "increasing") {
      message("increasing mode: retain marker with largest positive effect for heterozygous data")
      row_choice = ".N"
    } else if ((mode != "homozygous") & (mode != "increasing" | mode != "decreasing")) {
      stop("mode must be \n\t\"increasing\" or \"decreasing\" to use heterozygous biallelic data\n\t\"homozygous\" to keep only homozygous biallelic data.")
    }
    
    one = 1
    heterozygous = c("AGR", "CTY", "CGS", "ATW", "GTK", "ACM")
    
    if (mode == "homozygous") {
      effects = effects[num_occurrences == 2 & effect != 0]
    } else {
      number_of_markers = data.table::uniqueN(effects[num_occurrences == 3]$marker)
      progress <- utils::txtProgressBar(min = 0, max = number_of_markers, style = 3)
      effects = data.table::rbindlist(
        list(
          effects[num_occurrences == 2 & effect != 0],
          effects[num_occurrences == 3 , {
            alleles = sort(paste(.SD[, allele], collapse = ""))
            utils::setTxtProgressBar(progress, .GRP)
            if (alleles %in% heterozygous) {
              .SD[effect != 0][get(row_choice)]
            } else {
              NULL
            }
          }, by = marker][, c("marker", "trait", "marker_original", "chromosome", "position", "allele", "effect", "num_occurrences")]
        )
      )
    }    
    
    effects[, c("trait", "marker_original", "chromosome", "position", "allele", "num_occurrences") := NULL]
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
