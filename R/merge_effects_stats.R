# Setup and validate
#
# This function sets up and validates the data

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

merge_data <- function(stats_file, effects_file) {
  stats = read.table(stats_file, header=TRUE)
  effects = read.table(effects_file, header=TRUE)

  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>% group_by(Marker) %>% summarise(count=n()) %>% filter(count != 2)
  effects <- effects[!(effects$Marker %in% non_biallelic$Marker),]
  stats <- stats[!(stats$Marker %in% non_biallelic$Marker),]

  # Remove all NaN data to prevent math with NaN
  stats <- stats[!(stats$marker_F == "NaN"), ]

  # Split effects into even and odd rows and recombine into a single row without duplicate columns
  odd_effects<-effects[seq(1, nrow(effects), by = 2),]
  even_effects<-effects[seq(2, nrow(effects), by = 2),]
  effects <- left_join(odd_effects, even_effects, by = "Marker", "Trait")
  effects <- subset(effects, select = -c(Trait.y, Chr.y, Pos.y))

  # Merge stats and effects
  all_data <- left_join(stats, effects, by = "Marker")
  all_data <- subset(all_data, select = -c(Trait.x, Chr.x, Pos.x, dom_F, dom_p, add_F, add_p, marker_MS, error_MS, model_df, model_MS))

  # Delete temporary dataframes
  rm(non_biallelic)
  rm(even_effects)
  rm(odd_effects)

  # return merged data
  return(all_data)
}









