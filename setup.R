library(tidyverse)
setwd("/Users/Mason/Downloads")

effects = read.table("GrainColorGLMAllEffHiMAF.txt", header=TRUE)
stats = read.table("GrainColorGLMStatsHiMAF.txt", header=TRUE)

#Delete all markers in effects and stats with more or less alleles than 2
numbers <- effects %>% group_by(Marker) %>% summarise(count=n())
baddata <- numbers %>% filter(count != 2)
effects <- effects[!(effects$Marker %in% baddata$Marker),]
stats <- stats[!(stats$Marker %in% baddata$Marker),]

#Remove all NaN data due to it interfering with calculations
stats <- stats[!(stats$marker_F == "NaN"), ]

#Get all the even and odd number of rows in effects to later combine into one
odd_effects<-effects[seq(1, nrow(effects), by = 2),]
even_effects<-effects[seq(2, nrow(effects), by = 2),]

#Combine the odd and even rows of the effects datasheet and delete unecessary columns
CombinedEff <- left_join(odd_effects, even_effects, by = "Marker", "Trait")
CombinedEff <- subset(CombinedEff, select = -c(Trait.y, Chr.y, Pos.y))


#Combine Stats and CombinedEff to get all datasheets into one and delete unecessary columns
All <- left_join(stats, CombinedEff, by = "Marker")
All <- subset(All, select = -c(Trait.x, Chr.x, Pos.x))

#Going to go ahead and remove all unneeded datasheets
rm(numbers)
rm(even_effects)
rm(odd_effects)
rm(baddata)
rm(CombinedEff)