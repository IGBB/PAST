parse_SNP <- function(all_data, LD_file) {

  #Give upstream and downstream their own variables to make life easier
  LD_upstream <- LD_file[[1]]
  LD_downstream <- LD_file[[2]]

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
