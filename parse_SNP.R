library(tidyverse)
LD_file <- parse_LD("example/KernelColorLD.1.txt.gz")

#Give upstream and downstream their own variables to make life easier
LD_upstream <- LD_file[[1]]
LD_downstream <- LD_file[[2]]




for (name in names(LD_upstream)){
  temp_data <- LD_upstream[[name]]
  temp_data <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
}
LD_file <- LD_file[[1]][["1"]]
LD_file <- LD_file[(LD_file$Position1 %in% all_data$Pos),]
