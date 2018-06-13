library(tidyverse)
LD_file <- parse_LD("example/KernelColorLD.1.txt.gz")

LD_file1 <- LD_file[[1]][["1"]]
testfile <- LD_file[(LD_file[[1]][["1"]][["Position1"]] %in% all_data$Pos)]
