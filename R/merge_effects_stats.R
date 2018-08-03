merge_data <- function(stats_file, effects_file) {

  stats = read.table(stats_file, header = TRUE, sep = "\t")
  effects = read.table(effects_file, header = TRUE, sep = "\t")

  # Delete all markers in effects and stats with more or less alleles than 2
  non_biallelic <- effects %>% group_by(Marker) %>% dplyr::summarise(count=n()) %>% filter(count != 2)
  effects <- effects %>% filter(!(effects$Marker %in% non_biallelic$Marker), Locus != 0)
  stats <- stats %>% filter(!(stats$Marker %in% non_biallelic$Marker))

  # Remove all NaN data to prevent math with NaN
  stats <- stats %>% filter(MarkerR2 != "NaN", Chr != 0)

  # Split effects into even and odd rows and recombine into a single row without duplicate columns
  odd_effects<-effects[seq(1, nrow(effects), by = 2),]
  even_effects<-effects[seq(2, nrow(effects), by = 2),]
  effects <- merge(odd_effects, even_effects, by = "Marker") %>% mutate(Trait = Trait.x, Trait.x = NULL, Trait.y = NULL)

  # Delete temporary dataframes
  rm(non_biallelic)
  rm(even_effects)
  rm(odd_effects)

  # Merge stats and effects and return
  all_data <- merge(stats, effects, by = "Marker") %>% 
    mutate(Trait = Trait.x, Trait.x = NULL, Trait.y = NULL) %>% 
    select(-add_effect, -add_F, -add_p, -dom_F, -dom_p, -errordf)
}









