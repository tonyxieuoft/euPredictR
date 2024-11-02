
#' Parse Input Title
#'
#' Parse Input Description
#'
#' @param filename File name of the BLAST XML/json input file to parse
#' @param species Species from which the subject genome is from.
#'
#' @import rjson
#' @import dplyr
#' @import readr
#'
#' @export

parse_BLAST_json <- function(filename, species){

  # reads the json file, removes all newlines, then stores in list format
  json_data <- readr::read_file(filename) %>%
    gsub(pattern = "[\r\n]", replacement = "") %>%
    rjson::fromJSON()

  #instantiate a list to store all genes for the species
  genes <- list()

  for (query_report in json_data$BlastOutput2){

    query_title <- query_report$report$results$search$query_title

    # if there is no query title (incorrect file format)
    if (is.null(query_title)){
      stop("The input file is in the incorrect format")
    }
    else {
      # Continue
    }

    query_name <- strsplit(query_title, split = " ")[[1]][1]
    hits <- query_report$report$results$search$hits

    if (!is.list(hits)){
      stop("The input file is in the incorrect format")
    }
    else{
      # Continue
    }

    if (length(hits) < 1){
      genes[[query_name]] <- NULL
    }
    else{

      top_hit <- hits[[1]]
      hsps <- top_hit$hsps
      L <- length(hsps)

      hsp_data <- data.frame(q_start = numeric(length = L),
                             q_end = numeric(length = L),
                             q_seq = character(length = L),
                             s_start = numeric(length = L),
                             s_end = numeric(length = L),
                             s_seq = character(length = L),
                             s_strand = character(length = L),
                             e_val = numeric(length = L))

      for (hsp in top_hit$hsps){

        i <- hsp$num

        hsp_data$q_start[i] <- hsp$query_from
        hsp_data$q_end[i] <- hsp$query_to
        hsp_data$q_seq[i] <- hsp$qseq

        hsp_data$s_start[i] <- hsp$hit_from
        hsp_data$s_end[i] <- hsp$hit_to
        hsp_data$s_seq[i] <- hsp$hseq

        hsp_data$e_val[i] <- hsp$evalue
        hsp_data$s_strand[i] <- hsp$hit_strand
      }
      genes[[query_name]] <- hsp_data
    }
  }

  species_struct <- list(species = species, genes = genes)
  class(species_struct) <- "species"
  print(species_struct)
}
