#' Parse Linkage Disequilibrium
#'
#' @param LD_file The file containing linkage disequilibrium data (accepts plain-text or gzipped plain-text)
#' @importFrom rlang .data
#' @importFrom stats complete.cases
#' @import dplyr
#' @return The linkage disequilibrium data in a list containing dataframes for each chromosome.
#' @export
#'
#' @examples
#' demo_linkage_disequilibrium_file = system.file("extdata",
#' "LD.txt.xz", package = "PAST", mustWork = TRUE)
#' LD <- parse_LD(demo_linkage_disequilibrium_file)
parse_LD <- function(LD_file) {
  LD_all <- read.table(LD_file, header = TRUE) %>%
    dplyr::mutate(Dist_bp = ifelse(.data$Dist_bp == "N/A", NA, .data$Dist_bp))
  LD_all <-
    LD_all %>% dplyr::select(
      .data$Locus1,
      .data$Position1,
      .data$Site1,
      .data$Position2,
      .data$Site2,
      .data$Dist_bp,
      .data$R.2
    )
  LD_all <- LD_all[complete.cases(LD_all), ]
  split(LD_all, f = LD_all$Locus1)
}