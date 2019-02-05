parse_LD <- function(LD_file) {
  LD_all<-read.table(LD_file, header=TRUE) %>% mutate(Dist_bp = ifelse(Dist_bp == "N/A", NA, Dist_bp))

  LD_all<-LD_all %>% select("Locus1", "Position1", "Site1", "Position2", "Site2", "Dist_bp", "R.2")
  LD_all<-LD_all[complete.cases(LD_all),]

  # split by Locus1
  split(LD_all, f=LD_all$Locus1)
}