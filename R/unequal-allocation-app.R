#' R Shiny application to simulate clinical trials with restricted randomization 
#' @export

unequalAllocationApp <- function() {
  appDir <- system.file("shiny-apps", "unequal-allocation", package = "rartool")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `rartool`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
