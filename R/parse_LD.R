parse_LD <- function(LD_file) {
  LD_all<-read.table(LD_file, header=TRUE) %>% mutate(Dist_bp = ifelse(Dist_bp == "N/A", NA, Dist_bp))
  
  # filter out Locus1 == 0 and select certain columns and remove NaN
  # LD_all<-LD_all %>% filter(Locus1 == Locus2, Locus1 != 0) %>% 
  #   select("Locus1", "Position1", "Site1", "Position2", "Site2", "Dist_bp", "R.2")
  
  LD_all<-LD_all %>% select("Locus1", "Position1", "Site1", "Position2", "Site2", "Dist_bp", "R.2")
    
  LD_all<-LD_all[complete.cases(LD_all),]
  
  # split by Locus1
  LD<-split(LD_all, f=LD_all$Locus1)
  
  # already sorted by Site1
  LD_upstream<-LD
  
  # swap Site2/Position2 and Site1/Position1
  # for downstream, we treat SNP2 as SNP1
  LD_downstream<-LD
  for (name in names(LD_downstream)) {
    temp_data<-LD_downstream[[name]] %>% 
      mutate(temp_p2 = Position1, temp_s2 = Site1, Position1 = Position2, Site1 = Site2) %>%
      mutate(Site2 = temp_s2, Position2 = temp_p2, temp_s2 = NULL, temp_p2 = NULL)
    LD_downstream[[name]]<-temp_data[with(temp_data, order(Site1)),]
  }
  
  # return list with upstream and downstream
  list(LD_upstream, LD_downstream)
}