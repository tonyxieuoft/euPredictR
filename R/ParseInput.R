
#' Parse Input Title
#'
#' Parse Input Description
#'
#' @param filename file name of the BLAST XML/json input file to parse
#'
#' @import rjson
#' @import dplyr
#' @import readr
#'
#' @export



parse_json <- function(filename){
  list_data <- readr::read_file(filename) %>%
    gsub(pattern = "[\r\n]", replacement = "") %>%
    rjson::fromJSON()



  print(list_data$BlastOutput2)
}
