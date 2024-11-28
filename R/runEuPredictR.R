# Purpose: Run the shiny app for euPredictR
# Author: Tony Xie
# Date: 2024-11-27
# Version: 0.1.0
# Bugs and Issues: None

#' Launch Shiny App for euPredictR
#'
#' A function that launches the Shiny app for euPredictR. The app permits users
#' to upload directories containing BLAST output .json files, run the
#' package's prediction algorithm, and visualize the quality of the predictions
#' produced via a sequence coverage heatmap and gene-specific phylogenies.
#'
#' @return None, but opens up a shiny app
#'
#' @examples
#' \dontrun{
#' # Runs the shiny application for euPredictR
#' euPredictR::run_euPredictR()
#'
#' }
#'
#' @references
#' Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J,
#' McPherson J, Dipert A, Borges B. (2024). shiny: Web Application Framework
#' for R. R package version 1.9.1, https://CRAN.R-project.org/package=shiny.
#'
#' Attali, D., Sali, A. (2024). shinycssloaders: Add Loading Animations to a
#' 'shiny' Output While It's Recalculating. R package version 1.1.0,
#' https://CRAN.R-project.org/package=shinycssloaders.
#'
#' Pedersen, T., Nijs, V., Schaffner, T., Nantz, E. (2022). shinyFiles: A
#' Server-Side File System Viewer for Shiny. R package version 0.9.3,
#' https://CRAN.R-project.org/package=shinyFiles.
#'
#' @export
#' @importFrom shiny runApp
#' @importFrom shinycssloaders withSpinner
#' @import shinyFiles
#' @import shinyalert
#'
run_euPredictR <- function(){
  appDir <- system.file("shiny-scripts",
                        package = "euPredictR")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
