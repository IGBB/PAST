run_shiny <- function() {
  app_dir <- system.file("shiny-examples", "past", package = "past")
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}