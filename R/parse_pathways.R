parse_pathways <- function(LD, pathways_file) {
  
  # debugging code
  pathways_file = "example/pathways.txt"
  
  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")
  
  # merge pathways and LD by gene
  
  # filter out pathways with less than five genes associated with them
}