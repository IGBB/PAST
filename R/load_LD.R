#' Load Linkage Disequilibrium
#'
#' @param LD_file The file containing linkage disequilibrium data
#' @param LD_columns The names of the columns in your linkage disequilibrium
#'   data for the chromosome of the first SNP, the position of the first SNP,
#'   the site of the first SNP, the chromosome of the second SNP, the position
#'   of the second SNP, the site of the second SNP, the distance between the
#'   two SNPs, and the R.2
#' @importFrom rlang .data
#' @importFrom stats complete.cases
#' @import dplyr
#' @return The linkage disequilibrium data in a list containing
#'   dataframes for each chromosome.
#' @export
#'
#' @examples
#' demo_LD_file = system.file("extdata","LD.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' LD <- load_LD(demo_LD_file)
load_LD <- function(LD_file,
                    LD_columns = c("Locus1",
                                   "Position1",
                                   "Site1",
                                   "Position2",
                                   "Site2",
                                   "Dist_bp",
                                   "R.2")) {
  
  LD <- read.table(LD_file, header = TRUE, na.strings = c("N/A" , "NaN" )) %>%
    dplyr::mutate(Locus = as.character(!!as.name(LD_columns[1])),
                  Position1 = !!as.name(LD_columns[2]),
                  Site1 = !!as.name(LD_columns[3]),
                  Position2 = !!as.name(LD_columns[4]),
                  Site2 = !!as.name(LD_columns[5]),
                  Dist_bp = !!as.name(LD_columns[6]),
                  R.2 = as.numeric(!!as.name(LD_columns[7]))) %>% 
    dplyr::select(.data$Locus,
                  .data$Position1,
                  .data$Site1,
                  .data$Position2,
                  .data$Site2,
                  .data$Dist_bp,
                  .data$R.2)
  LD <- LD %>% filter(is.na(.data$R.2) != TRUE)
  split(LD, f = LD$Locus)
}