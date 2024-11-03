# and then also visualize output, no?


#' HeatMap Generator Helper
#'
#' Filler Description
#'
#' @param prediction_df Dataframe created by Build Predictions
get_seq_gaps <- function(prediction_df){

  gene_names <- get_unique_genes(prediction_df)
  species_names <- names(prediction_df)
  L <- length(gene_names)*length(species_names)

  seq_gaps_melted <- data.frame(species_name = character(length = L),
                                gene_names = character(length = L),
                                gap_percent = numeric(length = L))


  for (i in seq_along(prediction_df)){

    gene_names <- names(prediction_df[[i]])
    print(gene_names)
    for (seq in prediction_df[[i]]){
      print(seq)
    }
  }
}


#' Get coverage
#'
#' filler description
#'
#' @param seq The sequence to look for gaps for
#'
#' @import stringr
get_coverage <- function(seq){
  num_gaps <- stringr::str_count(seq, "-")
  seq_length <- nchar(seq)
  return((seq_length - num_gaps)/seq_length)
}


#' Get Unique Genes
#'
#' filler description
#'
#' @param prediction_df Database containing predictions, created by BuildPredictions
#'
get_unique_genes <- function(prediction_df){

  total_gene_list = new.env()
  for (genes in seq_along(prediction_df)){
    curr_gene_list <- names(genes)
    for (gene in curr_gene_list){
      if (is.null(total_gene_list[[gene]])){
        total_gene_list[[gene]] <- T
      }
    }
  }
  return(names(total_gene_list))
}
