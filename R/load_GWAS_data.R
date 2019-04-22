#' Load GWAS data
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param association_columns The names of the columns in your association data
#'   for Trait, Marker, Chromosome, Site, F, p, and marker_Rsquared
#' @param effects_columns The names of the columns in your effects data for
#'   Trait, Marker, Chromosome, Site, and effect
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @export
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_association_file = system.file("extdata", "association.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz",
#'   package = "PAST", mustWork = TRUE)
#' gwas_data <- load_GWAS_data(demo_association_file, demo_effects_file)
load_GWAS_data <- function(association_file,
                       effects_file,
                       association_columns = c("Trait",
                                         "Marker",
                                         "Locus",
                                         "Site",
                                         "F",
                                         "p",
                                         "marker_R2"),
                       effects_columns = c("Trait",
                                           "Marker",
                                           "Locus",
                                           "Site",
                                           "Effect")) {

    stats <- read.table(association_file, header = TRUE, sep = "\t") %>%
      dplyr::mutate(.data,
                  Trait = !!as.name(association_columns[1]),
                  Marker = !!as.name(association_columns[2]),
                  Chr = !!as.name(association_columns[3]),
                  Pos = !!as.name(association_columns[4]),
                  p = !!as.name(association_columns[5]),
                  marker_R2 = !!as.name(association_columns[6])) %>%
      dplyr::select(.data$Marker,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2)

    effects <- read.table(effects_file, header = TRUE, sep = "\t") %>%
      dplyr::mutate(.data,
                    Trait = !!as.name(effects_columns[1]),
                    Marker = !!as.name(effects_columns[2]),
                    Chr = !!as.name(effects_columns[3]),
                    Pos = !!as.name(effects_columns[4]),
                    Effect = !!as.name(effects_columns[5])) %>%
      dplyr::select(.data$Marker,
                    .data$Trait,
                    .data$Chr,
                    .data$Pos,
                    .data$Effect)

  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>%
    dplyr::group_by(.data$Marker) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  effects <-
    effects %>% dplyr::filter(!(.data$Marker %in% non_biallelic$Marker))
  stats <-
    stats %>% dplyr::filter(!(.data$Marker %in% non_biallelic$Marker))

  # Remove all NaN data to prevent math with NaN
  stats <- stats %>% dplyr::filter(.data$marker_R2 != "NaN")

  # Split effects into even and odd rows and
  # recombine into a single row without duplicate columns
  odd_effects <- effects[seq(1, nrow(effects), by = 2), ]
  even_effects <- effects[seq(2, nrow(effects), by = 2), ]
  effects <- merge(odd_effects, even_effects, by = "Marker")
  effects <- dplyr::mutate(
    effects,
    Trait = effects$Trait.x,
    Trait.x = NULL,
    Trait.y = NULL
  )

  # Merge stats and effects and return
  all_data <- merge(stats, effects, by = "Marker") %>%
    dplyr::mutate(
      Trait = .data$"Trait.x",
      Trait.x = NULL,
      Trait.y = NULL
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Chr,
      .data$Pos,
      .data$p,
      .data$marker_R2,
      .data$Effect.x,
      .data$Effect.y
    )
  all_data
}