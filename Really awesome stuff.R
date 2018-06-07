oddvals <- seq(1, nrow(GrainColorGLMAllEff), by = 2)
evenvals <- seq(2, nrow(GrainColorGLMAllEff), by = 2)
rows = nrow(GrainColorGLMAllEff)

for(i in 1:rows)
  merge(GrainColorGLMAllEff[oddvals[i], ], GrainColorGLMAllEff[evenvals[i], ], by = "Marker")


GrainColorGLMStats <- GrainColorGLMStats[!(GrainColorGLMStats$marker_F == "NaN"), ]
