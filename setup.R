library(tidyverse)
setwd("/Users/Mason/Downloads")

effects = read.table("GrainColorGLMAllEffHiMAF.txt", header=TRUE)
stats = read.table("GrainColorGLMStatsHiMAF.txt", header=TRUE)

#Get all the even and odd number of rows in effects to later combine into one
odd_effects<-effects[seq(1, nrow(effects), by = 2),]
even_effects<-effects[seq(2, nrow(effects), by = 2),]

#Combine the odd and even rows of the effects datasheet and delete unecessary columns
CombinedEff <- left_join(odd_effects, even_effects, by = "Marker", "Trait")
CombinedEff <- subset(CombinedEff, select = -c(Trait.y, Chr.y, Pos.y))


#Combine Stats and CombinedEff to get all datasheets into one and delete unecessary columns
All <- left_join(stats, CombinedEff, by = "Marker")
All <- subset(All, select = -c(Trait.x, Chr.x, Pos.x))

#Remove all NaN Data due to it interfering with Calculations
All <- All[!(All$marker_F == "NaN"), ]

#Going to go ahead and remove all unneeded datasheets
numbers <- effects %>% group_by(Marker) %>% summarise(count=n())

filteredstats <- stats[!(stats$marker_F == "NaN"), ]

effectstest <- effects[!()]