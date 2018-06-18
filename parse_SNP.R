library(tidyverse)

parse_SNP <- function(LD)
{
LD <- parse_LD("example/KernelColorLD.1.txt.gz")
all_data<-merge_data("../example/GrainColorGLMStatsHiMAF.txt", "../example/GrainColorGLMAllEffHiMAF.txt")

#Give upstream and downstream their own variables to make life easier
LD_upstream <- LD[[1]]
LD_downstream <- LD[[2]]



for (name in names(LD_upstream)){
  temp_data <- LD_upstream[[name]]
  LD_upstream[[name]] <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
}

for (name in names(LD_downstream)){
  temp_data <- LD_downstream[[name]]
  LD_downstream[[name]] <- temp_data[(temp_data$Position1 %in% all_data$Pos),]
}


list(LD_upstream, LD_downstream)
}

# all of the code below eventually needs to be moved to the for loops for upstream and downstream above
# for testing purposes, make sure you at least run the filtering loops above or the data below will be
# inaccurate

# get the data for the chromosome
chr1<-LD_upstream[["1"]]

# get the linked SNPs
chr1_linked<-LD_upstream[["1"]] %>% arrange(Position1) %>% filter(R.2 >= 0.8)

# get everything that didn't pass the filter
# only arrange by Position1 for downstream
# arrange by Position1 and Dist_bp for upstream
chr1_unlinked<-LD_upstream[["1"]] %>% filter(R.2 < 0.8) %>% arrange(Position1, Dist_bp)

# get the first instance of each unlinked SNP
chr1_unlinked<-chr1_unlinked[match(unique(chr1_unlinked$Position1), chr1_unlinked$Position1),]

# put linked and unlinked together and arrange by ascending Position1
chr1<-rbind(chr1_linked, chr1_unlinked) %>% arrange(Position1)

for (chr in names(LD_upstream)){
  new_data <- split(LD_upstream[[chr]], f=LD_upstream[[chr]]$Position1)
}