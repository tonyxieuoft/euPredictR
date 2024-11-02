# Test
#' Oxygen tag test (title)
#'
#' A test function that just prints out the name (description)
#'
#' @param name1 first item to be printed if we keep typing and typing and
#' @param name2 second item to be printed
#' @param name3 test filepath
#'
#' @return nothing!you can also see what specifically they should look like
#'
#' @export
#'

my_function <- function(name1, name2, name3){
  print(name1)
  print(name2)

  # remember in oxygen tags to also import necessities
  # namespace exports and imports automatically based on this
  #, which are generated through O2 tags
  # first line short title # remember, skipping line is v imp
  # second line longer explanation
  # and then just param
  # and then give some examples
  # to convert to rd. devtools::document()
  # you can define lists as a class if you want

  # add data to demonstrate how it works (goes in data directory)


  # usethis::use_data(Dataset) only after you convert it into a dataframe or
  # something else

  # for data tags, use same title description, source, describe
  # you need "Dataset" (of the name at the end)
  # what about fasta files...

  # 5 mb at most (probably just define a small dataset)

  # description, example, installation, overview

  # usethis::use_readme_rmd()
  # devtools::build_readme() -> just this once you've made edits
  # images must be in inst/extdata/...

  # read specifically needed sections in A4


  # remember just data.R and functions in R folder

  # use \dontrun{} if there's a problem (using dataset from another R package)

  # correct input -> correct output
  # easy test is if wrong input -> caught correctly

  # maybe just do an extremely simple test case


}
