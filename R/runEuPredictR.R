
#' Launch Shiny App for euPredictR
#'
#' A function that launches the Shiny app for euPredictR.
#'
#' @return No return value but open up a Shiny page
#'
#' @examples
#' \dontrun{
#' # Runs the shiny application for euPredictR
#' euPredictR::run_euPredictR()
#'
#' }
#'
#'
#' @export
#' @importFrom shiny runApp
#'
run_euPredictR <- function(){
  appDir <- system.file("shiny-scripts",
                        package = "euPredictR")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
