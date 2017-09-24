#' R Shiny application to simulate clinical trials with restricted randomization 
#' @export

restricted_app <- function() {
  appDir <- system.file("shiny-apps", "unequal-allocation", package = "artool")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `artool`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
