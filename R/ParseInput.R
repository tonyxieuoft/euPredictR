
#' Parse JSON-format BLAST output data
#'
#' This function facilitates the conversion of a BLAST output file in .json
#' format (in which the querie(s) are coding sequences and the set of subject
#' genomic sequences, ie. contigs/scaffolds/chromosomes, are from ONE species
#' only) into a nested R list data structure. The raw data in this format
#' serves as input for the build_predictions functions. This function assumes
#' the following:
#' \itemize{
#'  \item That the first word in the fasta heading(s) of the query coding
#'  sequence(s) is the abbreviated gene name for the query.
#'  \item That all exons for a gene are located on one genomic sequence, and not
#'  split between multiple contigs/scaffolds/chromosomes. Exons are generally
#'  close together and always on the same chromosome, so this should be
#'  reasonable.
#' }
#' @param filename File path of a single BLAST output file in .json format.
#' Output must be obtained in the context of input query coding sequences and a
#' set of genomic sequences from a SINGLE species.
#' @param species Name of species (string) represented by the set of subject
#' genomic sequences blasted against. Note that for the ease of the user, gene
#' names are extracted from query fasta headings rather than provided as a
#' direct argument
#' @param raw_blast_list S3 object of class RawBlastList, which is a nested R list
#' previously created by this parsing function. Contains raw BLAST output data.
#' Fill in only if you wish to group two sets of raw data together for
#' visualization.
#' The default value is NULL. See the results section below for more
#' information.
#'
#' @return Returns an S3 object of class RawBlastList, which stores the raw
#' BLAST output as a nested R list, as follows:
#' \itemize{
#'  \item The very outer list is in key-value format where each key is a
#'  species name and the corresponding value is a "gene list" for the species
#'  \item Each "gene list" is also in key-value format, where each key is a
#'  gene name and the  corresponding value is a list containing two attributes:
#'  1) \emph{hsp_data}, a data frame of high-scoring pairs identified by BLAST
#'  between the query coding sequence for the gene and subject genomic sequences
#'  and 2) \emph{query_len}, the length of the query coding sequence for the
#'  gene.
#'  \item Each row in \emph{hsp_data} corresponds to a single high-scoring pair,
#'  in which the columns are as follows:
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
#' If the parameter \emph{raw_blast_list} is not null, then the function adds
#' the newly curated data to raw_blast_list.
#'
#' @import dplyr
#'
#' @export
#'
#' @references
#'
#' Altschul S.F., Gish, W., Miller, W., Myers, E.W. and Lipman, D.J. (1990). Basic local alignment
#' search tool. Journal of Molecular Biology 215(3): 403-410.
#' https://doi.org/10.1016/S0022-2836(05)80360-2
#'
#' Wickham H, François R, Henry L, Müller K, Vaughan D (2023). _dplyr: A
#' Grammar of Data Manipulation_. R package version 1.1.4.
#' https://CRAN.R-project.org/package=dplyr.
#'
#' Wheeler, D. and Bhagwat, M. Blast QuickStart: Example-Driven Web-Based BLAST
#' Tutorial. In: Bergman NH, editor. Comparative Genomics: Volumes 1 and 2.
#' Totowa (NJ): Humana Press; 2007. Chapter 9.
#' https://www.ncbi.nlm.nih.gov/books/NBK1734/
#'
#' @examples
#'
#' # "orca_rho_BLAST_output.json" contains BLAST output for a rhodopsin query
#' # coding sequence and subject Orca genome
#' json_file1 <- system.file("extdata", "orca_rho_BLAST_output.json",
#'                            package="euPredictR")
#'
#' # "dolphin_rho_BLAST_output.json" contains BLAST output for a rhodopsin query
#' # coding sequence and subject dolphin (T. truncatus) genome
#' json_file2 <- system.file("extdata", "dolphin_rho_BLAST_output.json",
#'                            package="euPredictR")
#'
#'
#' # Example 1: Parse BLAST output data with no raw_blast_list input argument
#'
#' raw_blast_list1 <- parse_BLAST_json(filename = json_file1,
#'                                      species = "O_orca")
#'
#'
#' # Example 2: Parse BLAST output data and combine with raw blast list
#' # generated previously
#'
#' raw_blast_list2 <- parse_BLAST_json(filename = json_file2,
#'                                     species = "T_truncatus",
#'                                     raw_blast_list = raw_blast_list1)
#' # this list concatenates raw BLAST results for both orca and dolphin RHO

parse_BLAST_json <- function(filename, species, raw_blast_list=NULL){

  if (!is.character(species) | length(species) > 1){
    stop("Invalid input: species name must be a single string")
  }
  else {} # Continue

  if (!is.null(raw_blast_list) & class(raw_blast_list) != "RawBlastList"){
    stop("Invalid input: raw_blast_list must be in RawBlastList format and
         generated by a previous iteration of parse_BLAST_json.")
  }
  else{} # Continue

  # validate whether the file specified by 'filename' is in the correct format.
  # Also, convert the raw file into an R list with the exact same
  # attribute(key) to value structure. Uses a helper function as defiend below.
  json_data_as_list <- validate_BLAST_JSON_file(filename)

  #instantiate a list to store results for all genes for the species
  genes <- list()

  # Iterate through reports for the BLAST output that each contain results for
  # a single query. See validate_BLAST_JSON_file below for more detailed
  # documentation on the format.
  for (query_report in json_data_as_list$BlastOutput2){

    # get title of the query (validate_BLAST_JSON_file already checked to if
    # that this is possible)
    query_title <- query_report$report$results$search$query_title

    # get query length
    query_len <- query_report$report$results$search$query_len

    # take the first word in the query heading/title as the gene name
    # (assumption defined in description)
    gene_name <- strsplit(query_title, split = " ")[[1]][1]

    # get hits (each a collection of HSPs and represented by an HSP table) for
    # the query
    hits <- query_report$report$results$search$hits

    # if the query matches no hits, no HSP table is present to predict from for
    # the particular species-gene combination
    if (length(hits) < 1){
      genes[[gene_name]] <- NULL
    }
    else{

      # We assume that all exons for a gene are located on one genomic sequence.
      # Therefore, the top hit (with the best score) should contain HSPs
      # representing all of our exons of interest, and we discard the rest.
      top_hit <- hits[[1]]

      # Get high scoring pairs contained with the hit
      hsps <- top_hit$hsps

      # Convert the HSPs, which are in a nested list format, into a data frame
      # where each row is an HSP and the columns are HSP attributes (eg.
      # query/subject segments, boundaries, etc.). This is done by the following
      # steps:
      # 1) lapply converts each individual HSP into a data frame with one row
      # 2) do.call calls the rbind function to bind together all the individual
      # hsps into one data frame
      # 3) Remove the columns from the original dataset that we are not
      # interested in via select
      hsp_data <- do.call(rbind, lapply(hsps, data.frame)) %>%
        dplyr::select(-c(num, bit_score, score, identity, query_strand,
                           align_len, midline, gaps))

      # rename the columns of our HSP
      colnames(hsp_data) <- c("e_val", "q_start", "q_end", "s_start", "s_end",
                              "s_strand", "q_seq", "s_seq")

      # Add a seq_len column corresponding to the length of the HSP segment.
      # Then remove gaps in the HSP if the number in our query segment and
      # subject is equal (and therefore, only there for alignment purposes).
      hsp_data <- dplyr::mutate(hsp_data, seq_len = q_end - q_start + 1) %>%
        remove_gaps_if_equal() # helper function below

      # store the query length and HSP data from the top hit in our gene list
      genes[[gene_name]] <- list(query_len = query_len, hsp_data = hsp_data)
    }
  }

  # If raw_blast_list is a RawBlastList provided by the user, add the
  # species -> gene list key-value pair from the BLAST output .json file to it.
  if (!is.null(raw_blast_list)){
    raw_blast_list[[species]] <- genes
  }
  else{
    # otherwise, we create a new RawBlastList, and its first entry is the
    # species -> gene list key-value pair generated from the BLAST output .json
    # file
    raw_blast_list <- list()
    raw_blast_list[[species]] <- genes
    class(raw_blast_list) <- "RawBlastList"
  }
  return(raw_blast_list)
}


#' Validate BLAST JSON output file (Helper for parse_BLAST_json)
#'
#' Checks to see if the provided BLAST output file is 1) in .json format, and
#' 2) follows the conventions of results from NCBI's BLAST search tool. See the
#' example BLAST output files in the inst/extdata directory for more
#' information. In general, each BLAST output file must have the following
#' attribute structure:
#' \itemize{
#'  \item An overarching attribute 'BlastOutput2'. The value of 'BlastOutput2'
#'  is a list of 'reports', each containing results for a query.
#'  \item Within a 'report', the query title and length can be accessed by
#'  traversing the attributes: 'report' -> 'results' -> 'search'.
#'  Their respective attribute names are 'query_title' and 'query_len'
#'  \item The 'hits' attribute can also be accessed this way, and contains
#'  a list of hits each representing a collection of HSPs for a single subject
#'  sequence.
#' }
#'
#' @param filename File path to the BLAST JSON output file to validate.
#'
#' @return Return the contents of the BLAST JSON output file in R list format
#' as directly translated by the fromJSON() function from the rjson package.
#' The key-value 'dictionary' format of JSON is preserved in the list.
#'
#' @import configr
#' @import readr
#' @import rjson
#' @import dplyr
#'
#' @references
#'
#' Couture-Beil A (2024). rjson: JSON for R. R package
#' version 0.2.23.
#' https://CRAN.R-project.org/package=rjson.
#'
#' Li J. (2020). _configr: An Implementation of Parsing and
#' Writing Configuration File (JSON/INI/YAML/TOML)_. R
#' package version 0.3.5,
#' https://CRAN.R-project.org/package=configr.
#'
#' Wheeler, D. and Bhagwat, M. Blast QuickStart: Example-Driven Web-Based BLAST
#' Tutorial. In: Bergman NH, editor. Comparative Genomics: Volumes 1 and 2.
#' Totowa (NJ): Humana Press; 2007. Chapter 9.
#' https://www.ncbi.nlm.nih.gov/books/NBK1734/
#'
#' Wickham H, François R, Henry L, Müller K, Vaughan D. (2023). _dplyr: A
#' Grammar of Data Manipulation_. R package version 1.1.4.
#' https://CRAN.R-project.org/package=dplyr.
#'
#' Wickham H, Hester J, Bryan J. (2024). _readr: Read
#' Rectangular Text Data_. R package version 2.1.5.
#' https://CRAN.R-project.org/package=readr.
#'
#'
validate_BLAST_JSON_file <- function(filename){

  # First, check if file is in JSON format
  if (! configr::is.json.file(filename)) {
    stop("Inputted file is not in JSON format.")
  }
  else {} # Continue

  # reads the json file, removes all newlines, then stores in R list format
  json_data_as_list <- readr::read_file(filename) %>%
    gsub(pattern = "[\r\n]", replacement = "") %>%
    rjson::fromJSON()

  if (is.null(json_data_as_list$BlastOutput2)){
    stop("The input file is in the incorrect format: missing 'BlastOutput2'
         overarching attribute")
  }
  else{} # Continue

  # iterates through each report in the nested list to ensure correct structure
  # Note that each query_report corresponds to results for a single query
  for (query_report in json_data_as_list$BlastOutput2){

    # fasta heading of the query coding sequence, if format is correct
    query_title <- tryCatch(
      # wrapped in try catch because there are many points at which this could
      # fail, but all mean that the input is incorrect
      expr = {
        query_report$report$results$search$query_title
        },
      error = function(e){
        stop("The input file is in the incorrect format: missing BLAST
             attributes")
        }
    )

    # length of query coding sequence, if format is correct
    query_len <- tryCatch(
      # see reasoning above for try-catch implementation
      expr = {
        query_report$report$results$search$query_len
      },
      error = function(e){
        stop("The input file is in the incorrect format: missing BLAST
             attributes")
      }
    )

    # if there is no query title (incorrect file format)
    if (is.null(query_title)){
      stop("The input file is in the incorrect format: no query title")
    }
    else if (is.null(query_len)){
      stop("The input file is in the incorrect format: no query length")
    }
    else if (length(query_title) > 1 | !is.character(query_title)){
      stop("The input file is in the incorrect format:
           query header not a string")
    }
    else if (length(query_len) > 1 | !is.numeric(query_len)){
      stop("The input file is in the incorrect format: query length is not
           a number")
    }
    else{} # continue

    # note that each hit is the collection of hsps for a single subject sequence
    hits <- query_report$report$results$search$hits
    if (!is.list(hits)){
      stop("The input file is in the incorrect format: no hits list present")
    }
    else {} # Continue

    # Note that we do not check correctness of the HSP tables, as it is
    # inefficient to do so and it is highly unlikely at this point for user
    # input to be incorrect, regardless
  }
  return(json_data_as_list)
}


#' Remove Gaps in HSP Segment (Helper for parse_BLAST_json)
#'
#' Sometimes, high-scoring pairs (composed of a query segment and a subject
#' segment) contain gaps from alignment. When conducting predictions, it is
#' often useful to preserve these gaps to easily track 'missing' or 'extra'
#' nucleotides in the subject segment. However, if both the query and subject
#' have an equal number of gaps, removing them only shifts the alignment (there
#' aren't any "missing" or "extra" nucleotides), and helps remove any
#' frameshifting that might otherwise occur in the final gene prediction.
#'
#' @param hsp_data Data frame containing HSP data. Each row is an HSP, and the
#' columns are HSP attributes.
#'
#' @return Returns the altered HSP data frame, now with gaps removed if both
#' both the query and subject have an equal number of them.
#'
#' @import stringr
#'
#' @references
#'
#' Wheeler, D. and Bhagwat, M. Blast QuickStart: Example-Driven Web-Based BLAST
#' Tutorial. In: Bergman NH, editor. Comparative Genomics: Volumes 1 and 2.
#' Totowa (NJ): Humana Press; 2007. Chapter 9.
#' https://www.ncbi.nlm.nih.gov/books/NBK1734/
#'
#' Wickham H. (2023). stringr: Simple, Consistent Wrappers for Common String
#' Operations. R package version 1.5.1,
#' <https://CRAN.R-project.org/package=stringr>.
#'
remove_gaps_if_equal <- function(hsp_data){

  # iterates through each row of hsp_data
  for (i in seq_along(hsp_data$q_seq)){
    # count number of gaps in each segment
    gaps_q_seq <- stringr::str_count(hsp_data$q_seq[i], "-")
    gaps_s_seq <- stringr::str_count(hsp_data$s_seq[i], "-")

    # if the number of gaps between the two segments are equal, remove all gaps
    # from them
    if (gaps_q_seq == gaps_s_seq & gaps_q_seq > 0){
      hsp_data$q_seq[i] <- gsub('-', '', hsp_data$q_seq[i])
      hsp_data$s_seq[i] <- gsub('-', '', hsp_data$s_seq[i])
    }
  }
  return(hsp_data)
}

#' Parse Multiple BLAST JSON Output Files
#'
#' Given a directory containing multiple BLAST output files in .json format,
#' (where each BLAST run was conducted on coding sequences as queries and set
#' of subject genomic sequences for a species different from other runs),
#' generate a RawBlastList storing the data (see parse_BLAST_json for a more
#' detailed description). Assumptions are the same as parse_BLAST_json.
#'
#' @param dir_path Path to directory containing BLAST output .json files. Each
#' file must be the BLAST output from coding sequences as queries and
#' subject genomic sequences from a single species DIFFERENT from other runs.
#' The name of each .json file must be the species name for the subject set of
#' genomic sequences blasted against.
#' @param raw_blast_list Previous RawBlastList generated by the package. Use to
#' concatenate the newly generated RawBlastList and new one together. Default is
#' NULL.
#'
#' @return Returns an S3 object of type RawBlastList. See parse_BLAST_json for a
#' more detailed description.
#'
#' @import tools
#'
#' @export
#'
#' @references
#'
#' R Core Team. (2024). R: A Language and Environment for Statistical Computing.
#' R Foundation for Statistical Computing,
#' Vienna, Austria. https://www.R-project.org/.
#'
#' @examples
#'
#' # "sample_BLAST_output_dir" contains BLAST output for two different
#' # instances where a rhodopsin query is used against a subject set of Orca
#' # genomic sequences (results titled O_orca.json) and when it is used against
#' # a subject set of Dolphin genomic sequences (results titled
#' # T_truncatus.json)
#'
#' BLAST_out_directory <- system.dir("extdata", "example_raw_output_directory",
#'                                    package="euPredictR")
#'
#' # EXAMPLE: Convert data in BLAST_out_directory to RawBlastList format
#'
#' raw_blast_list1 <- parse_multiple_BLAST_json(dir_path = BLAST_out_directory,
#'                                              raw_blast_list = NULL)
#'
parse_multiple_BLAST_json <- function(dir_path, raw_blast_list=NULL){

  # get the list of files in the directory
  all_files <- list.files(path = dir_path, full.names = TRUE)

  for (filename in all_files){
    # get only the name of the file (no path), minus the extension
    species_name <- tools::file_path_sans_ext(basename(filename))
    # call parse_BLAST_json (essentially, parse_multiple_BLAST_json is just
    # a wrapper)
    raw_blast_list <- parse_BLAST_json(filename, species_name, raw_blast_list)
  }
  return(raw_blast_list)
}



# [END]
