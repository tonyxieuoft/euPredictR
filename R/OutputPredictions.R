
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
#' @returns Returns NULL. Output is either printed to the console or redirected
#' to the file path if provided.
#'
#' @examples
#'
#' # Get gene predictions from sample RawBlastList dataset provided in the
#' # package
#' gene_predictions <- build_predictions(sample_raw_blast_list)
#'
#' # Example 1: Output gene predictions in FASTA format to console as a string
#' \dontrun{
#'
#' output_predictions(predictions = gene_predictions,
#'                    output_file = NULL)
#' }
#'
#' # Example 2: Output gene predictions in FASTA format to a file
#' \dontrun{
#'
#' file_path <- /path/to/file # make sure there is no existing info in it
#' output_predictions(predictions = gene_predictions,
#'                    output_file = file_path)
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
output_predictions <- function(predictions, output_file = NULL){

  if (!inherits(predictions, "GenePredictions")){
    stop("Incorrect input. 'predictions' must be of type GenePredictions and
         produced by the build_predictions function in the package")
  }
  else{} # Continue

  if (!is.null(output_file) & !is.character(output_file)){
    stop("Incorrect input. Output file path must specified in string format.
         PLEASE NOTE that any information in the file (iif it exists) will be
         overwritten by this function.")
  }

  output_fasta <- ""

  # get all gene names in 'predictions' using the same methods as the
  # the heatmap creation function earlier
  total_gene_names <- reduce(lapply(predictions, names), dplyr::union)
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

# [END]
