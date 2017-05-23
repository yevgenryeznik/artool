unequalAllocationApp <- function() {
  appDir <- system.file("shiny-apps", "unequal-allocation", package = "rartool")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `rartool`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
