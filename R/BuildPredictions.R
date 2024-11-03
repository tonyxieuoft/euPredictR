
ALLOWED_OVERLAP <- 50
MAX_INTRON_LENGTH <- 1000000


#' Build Predictions with Preprocessed Data Frames
#'
#' Sample Description
#'
#' @param hsp_table Table of HSP values
#' @param query_len Length of the Query
#'
#' @import dplyr

build_prediction <- function(query_len, hsp_table){

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

#' Prediction wrapper
#'
#' sample description
#'
#' @param raw_carrier raw carrier preprocessed from before
#'
#' @export

build_predictions <- function(raw_carrier){

  #for (species in raw_carrier){
  #  for (gene in species){
  #    print(build_prediction(query_len = gene$query_len, hsp_table = gene$hsp_data))
  #  }
  #}
  prediction_carrier <- list()

  species_names <- names(raw_carrier)
  for (i in seq_along(raw_carrier)){
    # species -> looking at gene list

    new_gene_list <- list()

    gene_names <- names(raw_carrier[[i]])
    for (j in seq_along(raw_carrier[[i]])){
      gene <- raw_carrier[[i]][[j]]#$hsp_data
      predicted_seq <- build_prediction(query_len = gene$query_len,
                                        hsp_table = gene$hsp_data)

      new_gene_list[[gene_names[j]]] <- predicted_seq
    }

    prediction_carrier[[species_names[i]]] <- new_gene_list

  }
  return(prediction_carrier)
}

