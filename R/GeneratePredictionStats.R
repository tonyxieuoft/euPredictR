# Purpose: Generate and Plot Predicted Gene Coverage Data
# Author: Tony Xie
# Date: 2024-11-27
# Version: 0.1.0
# Bugs and Issues: None

#'
#' Create Gene Coverage Heatmap for Prediction Completeness
#'
#' Given an input nested list of gene predictions in GenePredictions format
#' (see documentation of the function 'build_predictions' for more), output
#' a heat map, where:
#' \itemize{
#' \item the x axis contains species names
#' \item the y axis contains gene names
#' \item the color gradient is the estimated coding sequence (CDS) coverage of
#' the predicted sequence for the given gene and species. Coding
#' sequence coverage is a measure of how complete the prediction is, and is
#' calculated by dividing the number of present nucleotides in the prediction
#' minus the expected length of the new coding sequence. The
#' more gaps present, the lower the coverage.
#' }
#'
#' @param predictions S3 object of class GenePredictions produced by an
#' iteration of the 'build_predictions' function in this package.
#'
#' @returns Returns a ggplot object representing the heatmap created.
#'
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
#' # make gene predictions from sample RawBlastList data provided in package
#' predictions <- build_predictions(sample_raw_blast_list,
#'                                  max_overlap = 30,
#'                                  max_intron_length = 50000)
#'
#' # Example: Create gene coverage/prediction completeness heatmap
#' \dontrun{
#'
#' gene_coverage_heatmap(predictions = predictions)
#'
#' }
#'
#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
#' New York, 2016.
#'
gene_coverage_heatmap <- function(predictions){

  if (!inherits(predictions, "GenePredictions")){
    stop("Incorrect input: 'predictions' must be in GenePredictions
         format and generated by the function build_predictions")
  }

  # calls helper function to get gene coverages in melted data frame format
  melted_coverages <- get_gene_coverage(predictions)

  # plot using ggplot
  coverage_heatmap <- ggplot2::ggplot(melted_coverages, aes(x = species_name,
                                       y = gene_name,
                                       fill = coverage,
                                       width = 0.9,
                                       height = 0.9)) +
    ggplot2::geom_tile() +
    # aesthetic options
    ggplot2::scale_fill_gradient(high = "blue4", low = "cornsilk3")+
    ggplot2::theme(panel.background = element_blank()) +
    ggplot2::xlab("Species Name") +
    ggplot2::ylab("Gene Name") +
    ggplot2::guides(fill=guide_legend(title="CDS coverage"),
                    x = guide_axis(angle = 270))

  return(coverage_heatmap)
}

#' Augment Dataset for Heatmap Generation (helper for gene_coverage_heatmap)
#'
#' Given an input GenePredictions object containing predicted coding sequences,
#' calculate gene coverage (ie. prediction completeness) for each predicted
#' sequence. Store the calculations in a melted data frame for heatmap
#' generation.
#'
#' @param predictions An S3 object of type GenePredictions, created by a
#' previous iteration of the build_predictions function in the package.
#'
#' @returns Return a data frame with three columns:
#' \itemize{
#'  \item species name
#'  \item gene name
#'  \item coding sequence coverage of the predicted coding sequence for the
#'  species and gene
#' }
#'
#' @import purrr
#'
#' @references
#' Wickham H, Henry L. (2023). purrr: Functional Programming Tools. R package
#' version 1.0.2. https://CRAN.R-project.org/package=purrr.
#'
get_gene_coverage <- function(predictions){

  # We get the total, unique gene names by the following:
  # 1) lapply applies the names function to each item in 'predictions'
  # to get the gene names present for each species
  # 2) reduce applies the union operator on the gene names until all
  # are concatenated into a single vector
  total_gene_names <- reduce(lapply(predictions, names), dplyr::union)

  # get species names
  species_names <- names(predictions)

  # Length of melted data frame to create is number of all species-genes
  # combinations
  L <- length(total_gene_names)*length(species_names)

  # instantiate melted data frame
  melted_coverages <- data.frame(species_name = character(length = L),
                                gene_name = character(length = L),
                                coverage = numeric(length = L))

  curr_index <- 1
  for (species in species_names){
    for (gene in total_gene_names){

      # get coding sequence prediction for a particular gene-species combination
      seq <- predictions[[species]][[gene]]

      melted_coverages$species_name[curr_index] <- species
      melted_coverages$gene_name[curr_index] <- gene

      # get coverage if sequence is present
      if (!is.null(seq)){
        melted_coverages$coverage[curr_index] <- calculate_coverage(seq)
      }
      else{
        # otherwise, coverage is assumed to be 0
        melted_coverages$coverage[curr_index] <- 0
      }
      curr_index = curr_index + 1
    }
  }

  return(melted_coverages)
}


#' Calculate Coding Sequence Coverage
#'
#' Calculate the percent of the subject species' coding sequence that was
#' successfully predicted for a particular gene. Does so by determining the
#' number of gaps in the sequence (corresponding to missing sections), then
#' dividing the remaining nucleotides by the total length.
#'
#' @param seq The predicted coding sequence to get sequence coverage for.
#'
#' @returns The coding sequence coverage of the predicted sequence, as a numeric
#' value
#'
#' @import stringr
#'
#' @references
#' Wickham H (2023). stringr: Simple, Consistent Wrappers for Common String
#' Operations. R package version 1.5.1.
#' https://CRAN.R-project.org/package=stringr.
#'
calculate_coverage <- function(seq){

  num_gaps <- stringr::str_count(seq, "-")
  seq_length <- nchar(seq)
  return((seq_length - num_gaps)/seq_length)
}


