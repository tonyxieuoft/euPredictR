
#' Parse JSON-format BLAST output data
#'
#' This function facilitates the conversion of a BLAST output file in .json
#' format (in which the quer(ies) are coding sequences and the set of subject
#' genomic sequences is from ONE species only) into a nested R list data
#' structure. The
#' raw data in this format serves as input for the build_predictions
#' functions. Assumes that the first word in the fasta heading(s) of the query
#' coding sequence(s) is the abbreviated gene name for the query.
#'
#' @param filename File path of a single BLAST output file in .json format.
#' Output must be obtained in the context of input query coding sequences and a
#' set of genomic sequences from a SINGLE species.
#' @param species Name of species (string) represented by the set of subject
#' genomic sequences blasted against. Note that for the ease of the user, gene
#' names are extracted from query fasta headings rather than provided as a direct argument
#' @param raw_list Nested R list previously created by the
#' package, contains raw BLAST output data. Fill in only if you wish to group
#' two sets of raw data together for visualization. THe default value is NULL.
#'
#' @return A nested R list data structure storing the raw BLAST output, as
#' follows:
#' \itemize{
#'  \item The outer list is in key-value format where each key is a
#'  species name and the corresponding value is a "gene list" for the species
#'  \item Each "gene list" is also in key-value format, where each key is a
#'  gene name and the  corresponding value is a list containing two attributes:
#'  1) \emph{hsp_data}, a data frame of high-scoring pairs identified by BLAST
#'  between the query coding sequence for the gene and subject genomic sequences
#'  and 2) \emph{query_len}, the length of the query coding sequence for the
#'  gene.
#'  \item Each row in \emph{hsp_data} corresponds to a single high-scoring pair,
#'  where the columns are as follows:
#'  \itemize{
#'    \item \emph{q_start}: start of the HSP segment on the query side.
#'    \item \emph{q_end}: end of the HSP segment on the query side.
#'    \item \emph{q_seq}: nucleotide sequence of the HSP segment on the query
#'    side
#'    \item \emph{s_start}, \emph{s_end} and \emph{s_seq} are the same, but on
#'    the subject side.
#'    \item \emph{s_strand}: The strand of the subject segment (Plus/Minus)
#'    \item \emph{e_val}: The HSP's e-value
#'    \item \emph{seq_len}: Length of the HSP
#'  }
#' }
#' The list is a concatenation of raw_list with the raw BLAST output pointed to
#' by filename.
#'
#'
#'
#' @import rjson
#' @import dplyr
#' @import readr
#' @import stringr
#' @import configr
#'
#' @export
#'
#' @examples
#' json_file1 <- system.file("extdata", "orca_rhfdsafo_BLAST_output.json",
#'                            package="euPredictR")
#' # Example 1: No raw_list input argument
#'
#' raw_list1 <- parse_BLAST_json(filename = json_file1, species = "O_orca")
#'
#' # Example 2: Using raw list generated previously as additional input
#'
#' raw_list2 <- parse_BLAST_json(filename = json_file1, species = "T_truncatus",
#'                               raw_list = raw_list1)
#' # this list concatenates raw BLAST results for both orca and dolphin RHO into
#' # one list

parse_BLAST_json <- function(filename, species, raw_list=NULL){

  if (! configr::is.json.file(filename)) {
    stop("Inputted file is not in JSON format.")
  }
  else {} # Continue

  # reads the json file, removes all newlines, then stores in list format
  json_data <- readr::read_file(filename) %>%
    gsub(pattern = "[\r\n]", replacement = "") %>%
    rjson::fromJSON()

  #instantiate a list to store all genes for the species
  genes <- list()

  for (query_report in json_data$BlastOutput2){

    query_title <- query_report$report$results$search$query_title
    query_len <- query_report$report$results$search$query_len

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

      hsp_data <- do.call(rbind, lapply(hsps, data.frame)) %>%
        subset(select = -c(num, bit_score, score, identity, query_strand,
                           align_len, midline, gaps))

      colnames(hsp_data) <- c("e_val", "q_start", "q_end", "s_start", "s_end",
                              "s_strand", "q_seq", "s_seq")

      hsp_data <- dplyr::mutate(hsp_data, seq_len = q_end - q_start + 1) %>%
        remove_gaps_if_equal()

      print(hsp_data)
      genes[[query_name]] <- list(query_len = query_len, hsp_data = hsp_data)
    }
  }

  if (!is.null(raw_list)){
    raw_list[[species]] <- genes
  }
  else{
    raw_list <- list()
    raw_list[[species]] <- genes
  }

  return(raw_list)
}

#' Title
#'
#' Description
validate_BLAST_JSON_file <- function(){
  print(1)
}



#' Parse Multiple BLAST json
#'
#' description
#'
#' @param dir_path path to directory containing json files from different
#' species
#' @param raw_list previous raw carrier, if available
#'
#' @import tools

parse_multiple_BLAST_json <- function(dir_path, raw_list=NULL){

  all_files <- list.files(path = dir_path, full.names = TRUE)

  raw_list <- NULL
  for (filename in all_files){
    species_name <- tools::file_path_sans_ext(basename(filename))
    raw_list <- parse_BLAST_json(filename, species_name, raw_list)
  }
  return(raw_list)
}


#' Check if number of gaps is equal and remove if so
#'
#' Dummy Description
#'
#' @param hsp_data Table containing HSP data
#'
#' @import stringr
#'
remove_gaps_if_equal <- function(hsp_data){

  for (i in seq_along(hsp_data$q_seq)){
    gaps_q_seq <- stringr::str_count(hsp_data$q_seq[i], "-")
    gaps_s_seq <- stringr::str_count(hsp_data$s_seq[i], "-")

    if (gaps_q_seq == gaps_s_seq & gaps_q_seq > 0){
      hsp_data$q_seq[i] <- gsub('-', '', hsp_data$q_seq[i])
      hsp_data$s_seq[i] <- gsub('-', '', hsp_data$s_seq[i])
    }
  }
  return(hsp_data)
}





