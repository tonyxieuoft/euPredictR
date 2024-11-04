# and then also visualize output, no?


#'
#' Create Heatmap
#'
#' filler description
#'
#' @param seq_gaps_melted Melted dataframe of sequence coverage
#'
#' @import ggplot2
create_seq_gap_heatmap <- function(seq_gaps_melted){
  ggplot2::ggplot(seq_gaps_melted, aes(x = species_name,
                                       y = gene_name,
                                       fill = coverage)) +
    geom_tile()
}


#' HeatMap Generator Helper
#'
#' Filler Description
#'
#' @param prediction_df Dataframe created by Build Predictions
get_seq_gaps <- function(prediction_df){

  total_gene_names <- get_unique_genes(prediction_df)
  print(total_gene_names)
  species_names <- names(prediction_df)
  L <- length(total_gene_names)*length(species_names)
  print(L)

  seq_gaps_melted <- data.frame(species_name = character(length = L),
                                gene_name = character(length = L),
                                coverage = numeric(length = L))

  curr_index <- 1

  for (i in seq_along(prediction_df)){

    gene_names <- names(prediction_df[[i]])

    for (j in seq_along(prediction_df[[i]])){
      seq_gaps_melted$species_name[curr_index] <- species_names[i]
      seq_gaps_melted$gene_name[curr_index] <- gene_names[j]
      seq_gaps_melted$coverage[curr_index] <-
        get_coverage(prediction_df[[i]][[j]])



      curr_index = curr_index + 1
    }

    non_included_genes <- dplyr::setdiff(gene_names, total_gene_names)
    for (gene in non_included_genes){

      seq_gaps_melted$species_name[curr_index] <- species_names[i]
      seq_gaps_melted$gene_name[curr_index] <- gene_names[j]
      seq_gaps_melted$coverage[curr_index] <- NA
      curr_index = curr_index + 1
    }
  }

  return(seq_gaps_melted)
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
  for (genes in prediction_df){
    curr_gene_list <- names(genes)
    for (gene in curr_gene_list){
      if (is.null(total_gene_list[[gene]])){
        total_gene_list[[gene]] <- T
      }
    }
  }
  return(names(total_gene_list))
}
