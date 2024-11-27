
#' Output BLAST predictions in FASTA format
#'
#' Given an input nested list of gene predictions in GenePrediction format,
#' output the predictions in fasta format. The output can be directly sent
#' to a file, or simply printed to the R console.
#'
#' @param predictions An S3 object of type GenePredictions, produced by the
#' package's build_predictions function.
#' @param output_file The file to output the fasta-format predictions to. NOTE
#' that anything in the output_file, if it exists will be OVERWRITTEN by this
#' function.
#'
#' @returns Returns NULL. Prints out the FASTA string representation of the
#' predicted sequence if no output file is provided; otherwise, redirects output
#' to the specified filepath.
#'
#' @examples
#'
#' # Get gene predictions from sample RawBlastList dataset provided in the
#' # package
#' gene_predictions <- build_predictions(sample_raw_blast_list,
#'                                       max_overlap = 30,
#'                                       max_intron_length = 50000)
#'
#' # Example 1: Output gene predictions in FASTA format to console as a string
#' \dontrun{
#'
#' output_predictions_fasta(predictions = gene_predictions,
#'                          output_file = NULL)
#' }
#'
#' # Example 2: Output gene predictions in FASTA format to a file
#' \dontrun{
#'
#' file_path <- /path/to/file # make sure there is no existing info in it
#' output_predictions_fasta(predictions = gene_predictions,
#'                          output_file = file_path)
#'
#' }
#' @export
#'
#' @import purrr
#'
#' @references
#' Wickham H, Henry L. (2023). purrr: Functional Programming Tools. R package
#' version 1.0.2. https://CRAN.R-project.org/package=purrr.
#'
output_predictions_fasta <- function(predictions, output_file = NULL){

  if (!inherits(predictions, "GenePredictions")){
    stop("Incorrect input. 'predictions' must be of type GenePredictions and
         produced by the build_predictions function in the package")
  }
  else{} # Continue

  if (!is.null(output_file) & !is.character(output_file)){
    stop("Incorrect input. Output file path must specified in string format.
         PLEASE NOTE that any information in the file (if it exists) will be
         overwritten by this function.")
  }

  output_fasta <- ""

  # get all gene names in 'predictions' using the same methods as the
  # the heatmap creation function earlier
  total_gene_names <- purrr::reduce(lapply(predictions, names), base::union)
  # get all species names
  species_names <- names(predictions)

  for (gene_name in total_gene_names){
    for (species in species_names){

      # get the gene prediction for the species/gene combination
      seq <- predictions[[species]][[gene_name]]

      # if it exists
      if (!is.null(seq)){

        # create a simple fasta heading composed of the species name and gene
        fasta_heading <- paste(">", gene_name, " ", species, "\n", sep = "")

        # add the heading and sequence to the growing fasta sequence
        output_fasta <- paste(output_fasta,
                              fasta_heading,
                              seq,
                              "\n\n",
                              sep = "")
      }
      else{} # Continue (don't print anything for this)
    }
  }

  # if no output file specified, print it out to the console
  if (is.null(output_file)){
    print(output_fasta)
  }
  else{
    # otherwise, write to specified file path
    write(output_fasta, output_file)
  }

  return(invisible(NULL))

}

#' Output predictions in data frame format
#'
#' Given an input nested list of gene predictions in GenePrediction format,
#' output the predictions as a data frame with three columns corresponding to
#' gene name, species and sequence.
#'
#' @param predictions  An S3 object of type GenePredictions, produced by the
#' package's build_predictions function.
#'
#' @returns Returns a data frame containing the predicted coding sequences. The
#' data frame has three string-containing columns: 'Gene', 'Species' and
#' 'Sequence', corresponding to the gene name, species and sequence for a given
#' prediction. If no coding sequence prediction is available for a gene-species
#' combination, no row for it is present.
#'
#' @examples
#'
#' # Get gene predictions from sample RawBlastList dataset provided in the
#' # package
#' gene_predictions <- build_predictions(sample_raw_blast_list,
#'                                       max_overlap = 30,
#'                                       max_intron_length = 50000)
#'
#' # Convert the output to a data frame containing the predicted coding
#' # sequences
#' output_predictions_data <- output_predictions_df(gene_predictions)
#'
#' @import purrr
#' @export
#'
#' @references
#' Wickham H, Henry L. (2023). purrr: Functional Programming Tools. R package
#' version 1.0.2. https://CRAN.R-project.org/package=purrr.
output_predictions_df <- function(predictions){

  if (!inherits(predictions, "GenePredictions")){
    stop("Incorrect input. 'predictions' must be of type GenePredictions and
         produced by the build_predictions function in the package")
  }

  # get all gene names and species names
  total_gene_names <- purrr::reduce(base::lapply(predictions, names),
                                    base::union)
  species_names <- names(predictions)

  # prepare an empty data frame to store all sequences
  num_sequences <- length(total_gene_names) * length(species_names)
  seq_dataframe <- data.frame("Gene" = character(length=num_sequences),
                              "Species" = character(length=num_sequences),
                              "Sequence" = character(length=num_sequences))
  j <- 1

  for (gene_name in total_gene_names){
    for (species in species_names){

      # get the gene prediction for the species/gene combination
      seq <- predictions[[species]][[gene_name]]

      # if it exists, populate the next row in the dataframe
      if (!is.null(seq)){

        seq_dataframe$Gene[j] <- gene_name
        seq_dataframe$Species[j] <- species
        seq_dataframe$Sequence[j] <- seq
        j = j+1
      }
    }
  }

  # remove empty columns and return dataframe
  return(seq_dataframe[seq_dataframe$Sequence != "",])
}

# [END]
