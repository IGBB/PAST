#' Load GWAS data
#'
#' @param files The file or files containing GWAS data as a vector
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @param input_type single, two, TASSEL
#' @param update_progress an optional function for use with shiny that updates the user on progress
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @export
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
#' @examples
#' demo_gwas_file = system.file("extdata", "single_gwas.txt.xz", package = "PAST", mustWork = TRUE)
#' gwas_data <- load_single_generic(c(demo_gwas_file))
#' demo_association_file = system.file("extdata", "association.txt.xz", package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.single_line.txt.xz", package = "PAST", mustWork = TRUE)
#' gwas_data <- load_two_generic(c(demo_association_file, demo_effects_file))
#' demo_association_file = system.file("extdata", "association.txt.xz", package = "PAST", mustWork = TRUE)
#' demo_effects_file = system.file("extdata", "effects.txt.xz", package = "PAST", mustWork = TRUE)
#' gwas_data <- load_TASSEL(c(demo_association_file, demo_effects_file))
load_GWAS_data <- function(files,
                           trait = "Trait",
                           marker = "Marker",
                           locus = "Locus",
                           site = "Site",
                           p = "p",
                           marker_R2 = "marker_R2",
                           effect = "Effect",
                           input_type,
                           effects_trait = "Trait",
                           effects_marker = "Marker",
                           effects_locus = "Locus",
                           effects_site = "Site",
                           update_progress = NULL) {
  
  if (length(files) == 1) {
    if (input_type == "single") {
      gwas_file = files[[1]]
      gwas_data = load_single_generic(gwas_file,
                                      trait,
                                      marker,
                                      locus,
                                      site,
                                      p,
                                      marker_R2,
                                      effect,
                                      update_progress)
    }
  } else if (length(files) == 2) {
    association_file = files[[1]]
    effects_file = files[[2]]
    if (input_type == "two") {
      gwas_data = load_two_generic(files[[1]], 
                                   files[[2]],
                                   trait,
                                   marker,
                                   locus,
                                   site,
                                   p,
                                   marker_R2,
                                   effect,
                                   effects_trait,
                                   effects_marker,
                                   effects_locus,
                                   effects_site,
                                   update_progress)
    } else if (input_type == "TASSEL") {
      gwas_data = load_TASSEL(files[[1]], 
                              files[[2]],
                              trait,
                              marker,
                              locus,
                              site,
                              p,
                              marker_R2,
                              effect,
                              effects_trait,
                              effects_marker,
                              effects_locus,
                              effects_site,
                              update_progress)
    }
  }
  gwas_data
}

#' Load single file data
#'
#' @param gwas_file The file containing GWAS data
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @param update_progress an optional function for use with shiny that updates the user on progress
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
load_single_generic <- function(gwas_file,
                                trait = "Trait",
                                marker = "Marker",
                                locus = "Locus",
                                site = "Site",
                                p = "p",
                                marker_R2 = "marker_R2",
                                effect = "Effect",
                                update_progress = NULL) {
  if (is.function(update_progress)) {
    parts = 2
  }
  
  if (is.function(update_progress)) {
    current_part = 0
    message = "Reading GWAS file"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  gwas_data <- read.table(gwas_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Marker = paste0(.data$Chr, "_", .data$Pos),
                  p = !!as.name(p),
                  marker_R2 = !!as.name(marker_R2),
                  Effect = !!as.name(effect)) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2,
                  .data$Effect)
  
  # Remove all NaN data to prevent math with NaN
  if (is.function(update_progress)) {
    current_part = 1
    message = "Removing NaN data"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  gwas_data <- gwas_data %>% dplyr::filter(.data$marker_R2 != "NaN")
  if (is.function(update_progress)) {
    current_part = 2
    message = "Complete"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  gwas_data %>% arrange(.data$Chr, .data$Pos)
}

#' Load two-file generic
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @param update_progress an optional function for use with shiny that updates the user on progress
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
load_two_generic <- function(association_file,
                             effects_file,
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
  
  if (is.function(update_progress)){
    parts = 4
  }
  
  if (is.function(update_progress)) {
    current_part = 0
    message = "Reading association file"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  stats <- read.table(association_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Marker = paste0(.data$Chr, "_", .data$Pos),
                  p = !!as.name(p),
                  marker_R2 = !!as.name(marker_R2)) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2)
  
  if (is.function(update_progress)) {
    current_part = 1
    message = "Reading effects file"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  effects <- read.table(effects_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(effects_trait),
                  Marker_original = !!as.name(effects_marker),
                  Chr = !!as.name(effects_locus),
                  Pos = !!as.name(effects_site),
                  Effect = !!as.name(effect)) %>%
    dplyr::select(.data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$Effect)
  
  # Remove all NaN data to prevent math with NaN
  if (is.function(update_progress)) {
    current_part = 2
    message = "Removing NaN data"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  stats <- stats %>% dplyr::filter(.data$marker_R2 != "NaN")
  
  # Merge stats and effects and return
  if (is.function(update_progress)) {
    current_part = 3
    message = "Combining association and effects data"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  all_data <- merge(stats, effects, by = "Marker_original") %>%
    dplyr::mutate(
      Trait = .data$"Trait.x",
      Trait.x = NULL,
      Trait.y = NULL
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Marker_original,
      Chr=.data$Chr.x,
      Pos=.data$Pos.x,
      .data$p,
      .data$marker_R2,
      .data$Effect
    )
  if (is.function(update_progress)) {
    current_part = 4
    message = "Complete"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  all_data %>% arrange(.data$Chr, .data$Pos)
}


#' Load TASSEL data
#'
#' @param association_file The association file
#' @param effects_file  The effects file
#' @param trait The names of the trait column
#' @param marker The names of the marker column
#' @param locus The names of the locus column
#' @param site The names of the site column
#' @param p The names of the p column
#' @param marker_R2 The names of the marker_R2 column
#' @param effect The names of the effect column
#' @param update_progress an optional function for use with shiny that updates the user on progress
#' @return The association data and the effects data merged into a dataframe
#'   with one row for each SNP
#' @importFrom rlang .data
#' @import utils
#' @import dplyr
load_TASSEL <- function(association_file,
                        effects_file,
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
  if (is.function(update_progress)) {
    parts = 5
  }
  
  if (is.function(update_progress)) {
    current_part = 0
    message = "Reading association file"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  stats <- read.table(association_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(trait),
                  Marker_original = !!as.name(marker),
                  Chr = !!as.name(locus),
                  Pos = !!as.name(site),
                  Marker = paste0(.data$Chr, "_", .data$Pos),
                  p = !!as.name(p),
                  marker_R2 = !!as.name(marker_R2)) %>%
    dplyr::select(.data$Marker,
                  .data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$p,
                  .data$marker_R2)
  
  if (is.function(update_progress)) {
    current_part = 1
    message = "Reading effects file"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  
  effects <- read.table(effects_file, header = TRUE, sep = "\t") %>%
    dplyr::mutate(Trait = !!as.name(effects_trait),
                  Marker_original = !!as.name(effects_marker),
                  Chr = !!as.name(effects_locus),
                  Pos = !!as.name(effects_site),
                  Effect = !!as.name(effect)) %>%
    dplyr::select(.data$Marker_original,
                  .data$Trait,
                  .data$Chr,
                  .data$Pos,
                  .data$Effect)
  
  # Delete all markers in effects and stats with more or less alleles than 2
  
  if (is.function(update_progress)) {
    current_part = 2
    message = "Removing non-biallelic data"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  non_biallelic <- effects %>%
    dplyr::group_by(.data$Marker_original) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::filter(count != 2)
  effects <-
    effects %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  stats <-
    stats %>% 
    dplyr::filter(!(.data$Marker_original %in% non_biallelic$Marker_original))
  
  
  # Remove all NaN data to prevent math with NaN
  if (is.function(update_progress)) {
    current_part = 3
    message = "Removing NaN data"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  stats <- stats %>% dplyr::filter(.data$marker_R2 != "NaN")
  
  # Split effects into even and odd rows and
  # recombine into a single row without duplicate columns
  if (is.function(update_progress)) {
    parts = 5
    current_part = 4
    message = "Combining association and effects data"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  odd_effects <- effects[seq(1, nrow(effects), by = 2), ]
  even_effects <- effects[seq(2, nrow(effects), by = 2), ]
  effects <- merge(odd_effects, even_effects, by = "Marker_original")
  effects <- dplyr::mutate(
    effects,
    Trait = effects$Trait.x,
    Trait.x = NULL,
    Trait.y = NULL
  )
  
  # Merge stats and effects and return
  all_data <- merge(stats, effects, by = "Marker_original") %>%
    dplyr::mutate(
      Trait = .data$"Trait.x",
      Trait.x = NULL,
      Trait.y = NULL
    ) %>%
    dplyr::select(
      .data$Marker,
      .data$Marker_original,
      .data$Chr,
      .data$Pos,
      .data$p,
      .data$marker_R2,
      Effect=.data$Effect.x
    )
  
  if (is.function(update_progress)) {
    current_part = 5
    message = "Complete"
    update_progress(message = message, value = 100/parts/100*current_part, paste0(round(100/parts*current_part, 2), "%"))
  }
  
  all_data %>% arrange(.data$Chr, .data$Pos)
}
