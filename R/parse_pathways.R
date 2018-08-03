parse_pathways <- function(tagSNPs, pathways_file) {
  
  # debugging code
  pathways_file = "example/pathways.txt"
  
  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")
  
  # UP/DOWNSTREAM LOOP
  for (i in 1:length(tagSNPs)) {
    tagSNPs_stream <- tagSNPs[[i]]
    names(LD_stream)
    
    # BEGIN PROCESSING BY CHROMOSOMES LOOP
    for (name in names(tagSNPs_stream)) {
      temp_tagSNPs <-tagSNPs_stream[[name]] 
      
      # merge pathways and temp_tagSNPs by gene
      # important functions: merge
      # reference: parse_SNP.R, line 147
      
      # filter out pathways with less than five genes associated with them
      # important functions: group_by, summary, filter
      # reference: parse_SNP.R, line 118-120
      
    } # END CHROMOSOMES LOOP
    
  } # END UP/DOWNSTREAM LOOP
  
  
  
  
}