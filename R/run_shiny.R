#' Run Shiny App
#'
#' @import dplyr
#' @import ggplot2
#' @import shiny
#' @import shinydashboard
#' @importFrom gridExtra grid.arrange
#' @importFrom DT renderDataTable
#' @importFrom DT dataTableOutput
#' @export
#' @return Does not return a value
run_shiny <- function() {
  app_dir <- system.file("shiny-examples", "past", package = "PAST")
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `PAST`.",
         call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}
