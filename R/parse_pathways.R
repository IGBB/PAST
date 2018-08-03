parse_pathways <- function(tagSNPs, pathways_file) {
  
  # debugging code to be removed
  pathways_file = "example/pathways.txt"
  
  # load pathways
  pathways <- read.table(pathways_file, sep = "\t", header = TRUE, quote="")
  
  # UP/DOWNSTREAM LOOP
  for (i in 1:length(tagSNPs)) {
    tagSNPs_stream <- tagSNPs[[i]]

    # BEGIN PROCESSING BY CHROMOSOMES LOOP
    for (name in names(tagSNPs_stream)) {
      temp_tagSNPs <-tagSNPs_stream[[name]] 
      
      # merge pathways and temp_tagSNPs by gene
      # important functions: merge
      # reference: parse_SNP.R, line 147
      
      # filter out pathways with less than five genes associated with them
      # important functions: group_by, summary, filter
      # reference: parse_SNP.R, line 118-120
      
      # store pathway information for chromosomes
      tagSNPs_stream[[name]] <- temp_tagSNPs
      
    } # END CHROMOSOMES LOOP
    
    # store pathway information for up/dowstream
    tagSNPs[[i]] <- tagSNPs_stream
    
  } # END UP/DOWNSTREAM LOOP
  
# return updated tagSNPs with pathway information
tagSNPs

}