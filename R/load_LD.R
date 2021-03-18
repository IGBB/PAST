#' Load Linkage Disequilibrium
#'
#' @param LD_file the file containing linkage disequilibrium data
#' @param locus the name of the column containing the locus or chromosome
#' @param position1 the name of the column containing the first SNP's position
#' @param site1 the name of the column containing the first SNP's site
#' @param position2 the name of the column containing second SNP's position
#' @param site2 the name of the column containing the second SNP's site
#' @param distance the name of the column containing distance between the SNPs
#' @param r_squared the name of the column containing R^2
#' @return the linkage disequilibrium data in a data.table
#' @export
#' @import data.table
#' @examples
#' LD_file = system.file("extdata","LD.txt.gz", package = "PAST", 
#' mustWork = TRUE)
#' LD <- load_LD(LD_file)
load_LD <- function(LD_file,
                    locus = "Locus1",
                    position1 = "Position1",
                    site1 = "Site1",
                    position2 = "Position2",
                    site2 = "Site2",
                    distance = "Dist_bp",
                    r_squared = "R^2") {
  LD_columns = c(locus, position1, site1, position2, site2, distance, r_squared)
  
  # Read the file and select its columns.
  # Set the names to values used throughout PAST instead of what the user
  #   provided.
  # Make the locus column a character, just in case it isn't.
  # Drop rows with NA values in the distance or r_squared columns.
  LD <- data.table::fread(LD_file, 
                          select = LD_columns, 
                          na.strings = c("N/A" , "NaN" ))
  data.table::setnames(LD, LD_columns, c("locus",
                                         "position1",
                                         "site1",
                                         "position2",
                                         "site2",
                                         "distance",
                                         "r_squared"))
  LD[, locus:=as.character(locus)]
  LD <- na.omit(LD, cols = c("distance", "r_squared"))
  return(LD)
}