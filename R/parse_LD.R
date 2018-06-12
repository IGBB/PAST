# Setup and validate
#
# This function sets up and validates the data

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

parse_LD <- function(LD_file) {
  LD_all<-read.table(LD_file, header=TRUE)

  # filter out Locus1 == 0 and select certain columns and remove NaN
  LD_all<-filter(LD_all, Locus1 != 0) %>% select("Locus1", "Position1", "Site1", "Position2", "Site2", "Dist_bp", "R.2")
  LD_all<-LD_all[complete.cases(LD_all),]

  # split by Locus1
  LD<-split(LD_all, f=LD_all$Locus1)

  # already sorted by Site1
  LD_upstream<-LD

  # sort by Site2 for downstream
  LD_downstream<-LD
  for (name in names(LD_downstream)) {
    temp_data<-LD_downstream[[name]]
    LD_downstream[[name]]<-temp_data[with(temp_data, order(Site2)),]
  }

  # return list with upstream and downstream
  list(LD_upstream, LD_downstream)
}









