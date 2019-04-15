#' Merge Data
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param type The type of analysis from TASSEL (mlm or glm)
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @export
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_association_file = system.file("extdata",
#' "association.txt.xz", package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata",
#' "effects.txt.xz", package = "PAST", mustWork = TRUE)
#' type = "mlm"
#' merged_data <- merge_data(demo_association_file, demo_effects_file, type)
merge_data <- function(association_file, effects_file, type) {
  stats <- read.table(association_file, header = TRUE, sep = "\t")
  effects <- read.table(effects_file, header = TRUE, sep = "\t")

  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>%
    dplyr::group_by(.data$Marker) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  effects <-
    effects %>% dplyr::filter(!(.data$Marker %in% non_biallelic$Marker))
  stats <-
    stats %>% dplyr::filter(!(.data$Marker %in% non_biallelic$Marker))

  if (type == "mlm") {
    stats <- stats %>% dplyr::mutate(Chr = .data$Locus, Pos = .data$Site)
  } else if (type == "glm") {
    stats <- stats %>% dplyr::mutate(marker_R2 = .data$marker_Rsq, F = .data$marker_F)
    effects <- effects %>% dplyr::mutate(Effect = .data$Estimate)
  } else if (type == "other"){
    #
  } else {
    print("Type not specified correctly. Type must be \"mlm\" or \"glm\".")
    return()
  }

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
      Trait.y = NULL,
      Marker = paste0("S", .data$Chr, "_", .data$Pos)
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Chr,
      .data$Pos,
      .data$F,
      .data$p,
      .data$marker_R2,
      # .data$Allele.x,
      .data$Effect.x,
      # .data$Obs.x,
      # .data$Allele.y,
      .data$Effect.y,
      # .data$Obs.y,
      .data$Trait
    )
  all_data
}