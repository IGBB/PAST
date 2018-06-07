# Include libraries so that others know what dependencies to load
library(tidyverse)

# Include data loading as generic names, since not every experiement will be grain color
# We'll change the way this works later to be user input
effects = read.table("GrainColorGLMAllEffHiMAF.txt", header=TRUE)
stats = read.table("GrainColorGLMStatsHiMAF.txt", header=TRUE)


# This array has all the odd and even integers. You can reduce some of your code by using it like I did below.
oddvals <- seq(1, nrow(effects), by = 2)
evenvals <- seq(2, nrow(effects), by = 2)
rows = nrow(effects)

# In general, you can do most things in R without loops. You can get the odd values by doing what I've written below.
for(i in 1:rows)
  merge(effects[oddvals[i], ], effects[evenvals[i], ], by = "Marker")

# Get odd lines from effects and even lines from effects
odd_effects<-effects[seq(1, nrow(effects), by = 2),]
even_effects<-effects[seq(2, nrow(effects), by = 2),]

GrainColorGLMStats <- GrainColorGLMStats[!(GrainColorGLMStats$marker_F == "NaN"), ]
