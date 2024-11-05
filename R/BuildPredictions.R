
#' Build Gene Predictions from BLAST output data.
#'
#' Given that BLAST output data has been pre-processed by the parse_BLAST_json
#' or parse_multiple_BLAST_json functions, this function takes in data in
#' RawBlastList format, predicts protein-coding sequences via dynamic
#' programming algorithms and outputs their sequences in a nested list
#' (elaborated further in the results section). On a high level, for each gene,
#' the prediction algorithm identifies compatible HSPs with the greatest
#' gene coverage, sorts out overlaps between them, and finally stitches them
#' together to create a new coding sequence (complete or incomplete). Note that
#' previous assumptions stated for the functions that parse input hold here as
#' well (the one that states all relevant HSPs are on the same genomic subject
#' sequence is especially relevant).
#' HSP compatibility (whether two HSPs can fit into the same gene
#' model) is based on the additional assumption that relevant HSPs each loosely
#' correspond to an exon in the gene. From this, we draw the following
#' observations for compatible HSPs 'A' and 'B':
#' \itemize{
#'  \item \emph{Location Linearity}: If the subject segment of 'A' is
#'  before 'B' in the genomic sequence, the query segment of 'A' must
#'  also be before the query segment of 'B'.
#'  \item \emph{Maximum Overlap Restriction}: The query segment for 'A' must
#'  overlap minimally with the query segment for 'B'(a higher threshold above
#'  0 is set for some wiggle room).
#'  \item \emph{Maximum Intron Length}: The subject segment for 'A' must be
#'  somewhat close in proximity to 'B', otherwise the presumed intron(s) between
#'  them would be massively long.
#'  \item \emph{Same Strand}: The subject segment for 'A' must be on the strand
#'  as the subject segment for 'B' (plus/minus).
#' }
#' The maximum intron length and maximum overlap can be modified in the
#' parameter settings of the function.
#'
#' @param raw_blast_list S3 object of type RawBlastList, containing BLAST output
#' data generated from query coding sequences and subject genomic sequences.
#' @param max_overlap A non-negative integer representing the maximum overlap
#' between the query segments of two compatible HSPs. The default is set to 30.
#' @param max_intron_length A non-negative integer representing the maximum
#' distance between the subject segments of two adjacent, compatible HSPs. The
#' default is set to 50000.
#'
#' @returns Returns an S3 object of type GenePredictions, which is a
#' nested list containing predicted protein-coding sequences with the following
#' structure:
#' \itemize{
#' \item The outer list follows the key-value format,
#' where species names are the keys and each value is a 'gene list' for the
#' species. (just like RawBlastList)
#' \item Each 'gene list' is also a list in key-value format, where each key is
#' a gene name and each value is a predicted sequence in string format. Missing
#' sections of the predicted coding sequence are denoted "-".
#' }
#'
#' @examples
#' # example code
#'
#'
#' @export
build_predictions <- function(raw_blast_list,
                              max_overlap=30,
                              max_intron_length=50000){

  # instantiate a list that will be defined as a GenePredictions object later
  predictions <- list()

  # names of represented species to predict sequences for
  species_names <- names(raw_blast_list)

  # iterates along indices corresponding to items where the key is a species
  # name and value is a gene list
  for (i in seq_along(raw_blast_list)){

    # instantiate a new gene list that will hold gene name -> sequence
    # key-value pairs as items
    new_gene_list <- list()

    # list of genes to predict for a particular species
    gene_names <- names(raw_blast_list[[i]])

    # iterates along indices corresponding to items where the key is a gene name
    # and value is a list containing a query's length and its HSP table
    for (j in seq_along(raw_blast_list[[i]])){

      gene <- raw_blast_list[[i]][[j]] # the hsp_table/query_len list

      # Predict a single coding sequence based on the HSP table, for the given
      # gene and species
      predicted_seq <- build_prediction(query_len = gene$query_len,
                                        hsp_table = gene$hsp_data,
                                        max_overlap = max_overlap,
                                        max_intron_length = max_intron_length)

      # add the sequence to our growing gene prediction list (for one species)
      new_gene_list[[gene_names[j]]] <- predicted_seq
    }

    # add the new gene list containing predicted sequences for one species to
    # the overarching predictions list
    predictions[[species_names[i]]] <- new_gene_list

  }
  class(predictions) <- "GenePredictions"
  return(predictions)
}


#' Build Predictions with Preprocessed Data Frames
#'
#' Sample Description
#'
#' @param hsp_table Table of HSP values
#' @param query_len Length of the Query
#'
#' @import dplyr

build_prediction <- function(query_len,
                             hsp_table,
                             max_overlap,
                             max_intron_length){

  if (is.null(hsp_table)){
    return(NULL)
  }
  else{
    # Continue
  }

  sorted_table <- dplyr::arrange(hsp_table, q_start)
  num_hsps <- nrow(sorted_table)

  dp <- data.frame(parent = numeric(length = num_hsps),
                   coverage = numeric(length = num_hsps))

  for (i in 1:num_hsps){

    dp$parent[i] <- -1
    dp$coverage[i] <- sorted_table$seq_len[i]

    max_coverage <- dp$coverage[i]

    for (j in 1:i){
      if (is_compatibile(sorted_table, j, i)){
        overlap_with_j <- max(0,
                              sorted_table$q_end[j] - sorted_table$q_start[i]+1)
        new_coverage <- dp$coverage[j] +
          sorted_table$seq_len[i] -
          overlap_with_j
        if (new_coverage > max_coverage){
          max_coverage <- new_coverage
          dp$parent[i] <- j
        }
      }
    }
    dp$coverage[i] <- max_coverage
  }

  #print(sorted_table)
  #print(dp)
  return(create_sequence(dp, sorted_table, query_len))
}


#' Create Sequence from DP table
#'
#' Dummy Description
#'
#' @param dp The dynamic programming table
#' @param hsp_table The table of hsps
#' @param query_len Length of the query
create_sequence <- function(dp, hsp_table, query_len){

  ending_index <- which.max(dp$coverage)
  seq <- recursive_create(dp, hsp_table, ending_index)

  ending_gap <- strrep("-", query_len - hsp_table$q_end[ending_index])
  return(paste(seq, ending_gap, sep = ""))
}

#' Recursive Helper
#'
#' Dummy Description
#'
#' @param dp The dynamic programming table
#' @param hsps The table of hsps
#' @param i index of the last HSP in the sequence to create
#'
#' @import stringr
recursive_create <- function(dp, hsps, i){

  if (dp$parent[i] == -1){
    missing_beginning <- strrep("-", hsps$q_start[i] - 1)
    return(paste(missing_beginning, hsps$s_seq[i], sep = ""))
  }
  else{

    p <- dp$parent[i]
    seq <- recursive_create(dp, hsps, p)

    if (hsps$q_end[p] < hsps$q_start[i]){

      gap <- strrep("-", hsps$q_start[i] - hsps$q_end[p] - 1)
      return(paste(seq, gap, hsps$s_seq[i], sep = ""))
    }
    else{
      return(handle_overlap(hsps, seq, p, i))
    }
  }
}

#' Merge Overlap Helper
#'
#' Placeholder description
#'
#' @param hsps Table of HSPs
#' @param seq Currently built sequence
#' @param p Index of parent sequence
#' @param i Index of current sequence
#'
handle_overlap <- function(hsps, seq, p, i){

  overlap_num <- hsps$q_end[p] - hsps$q_start[i] + 1

  num_left <- overlap_num
  p_seq <- hsps$q_seq[p]
  j <- nchar(p_seq)

  while (num_left > 0){
    if (substr(p_seq, j, j) != "-"){
      num_left <- num_left - 1
    }
    else{
      # continue
    }
    j <- j - 1
  }

  seq_low_bound <- nchar(seq) - (nchar(p_seq) - (j + 1))
  seq_overlap <- substr(seq, seq_low_bound, nchar(seq))

  num_right <- overlap_num
  i_seq <- hsps$q_seq[i]
  curr <- 1

  while (num_right > 0){
    if (substr(i_seq, curr, curr) != "-"){
      num_right <- num_right - 1
    }
    else{
      # continue
    }
    curr <- curr + 1
  }

  i_bound <- curr - 1
  i_overlap <- substr(hsps$s_seq[i], 1, i_bound)


  if (stringr::str_count(seq_overlap, "-") == 0){
    i_to_take <- substr(hsps$s_seq[i], i_bound + 1, nchar(hsps$s_seq[i]))
    return(paste(seq, i_to_take, sep = ""))
  }
  else{
    seq_to_take <- substr(seq, 1, seq_low_bound-1)
    return(paste(seq_to_take, hsps$s_seq[i], sep = ""))
  }
}


#' Helper function for building predictions.
#'
#' Determines if two HSPs are compatible
#'
#' @param hsp_table a sorted table of HSPs
#' @param j index of HSP with lower query start
#' @param i index of HSP with high query start
#'
is_compatibile <- function(hsp_table, j, i){

  # if the query seq for j overlaps with all of i
  if (hsp_table$q_end[j] >= hsp_table$q_end[i]){
    return(F)
  }
  # if query seq for j overlaps with i by more than ALLOWED_OVERLAP base pairs
  else if (hsp_table$q_end[j] - hsp_table$q_start[i] > ALLOWED_OVERLAP){
    return(F)
  }
  # if the subject sequences for i and j are not on the same strand
  else if (hsp_table$s_strand[j] != hsp_table$s_strand[i]){
    return(F)
  }
  # check for hsp linearity
  else if (hsp_table$s_start[j] < hsp_table$s_end[j] &
      hsp_table$s_end[j] >= hsp_table$s_start[i]){
    return(F)
  }
  # check for hsp linearity
  else if (hsp_table$s_start[j] > hsp_table$s_end[j] &
           hsp_table$s_end[j] <= hsp_table$s_start[i]){
    return(F)
  }

  else if (abs(hsp_table$s_end[j] - hsp_table$s_start[i]) > MAX_INTRON_LENGTH){
    return(F)
  }
  else{
    return(T)
  }
}

