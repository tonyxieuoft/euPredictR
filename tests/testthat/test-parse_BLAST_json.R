#context("Test the parse_BLAST_json file in ParseInput.R")
library(euPredictR)


test_that("Valid .json file for dolphin and orca IL6,IL10,TNF,RHO, SOX2 HSPs
          obtained by
          BLAST (see data.R for details)", {
   good_filename1 <- system.file("extdata/example_raw_output_directory",
                                "T_truncatus.json",
                                package = "euPredictR")
   good_filename2 <- system.file("extdata/example_raw_output_directory",
                                 "O_orca.json",
                                 package = "euPredictR")

   raw_blast_list1 <- parse_BLAST_json(
     filename = good_filename1,
     species = "T_truncatus",
     raw_blast_list = NULL
   )

   raw_blast_list2 <- parse_BLAST_json(
     filename = good_filename2,
     species = "O_truncatus",
     raw_blast_list = raw_blast_list1
   )

   testthat::expect_s3_class(raw_blast_list1, "RawBlastList")
   testthat::expect_s3_class(raw_blast_list2, "RawBlastList")
})


test_that("Checking for invalid user input for parse_BLAST_json", {

  good_filename <- system.file("extdata/example_raw_output_directory",
                               "T_truncatus.json",
                               package = "euPredictR")
  testthat::expect_error(parsed_output <- parse_BLAST_json(
    filename = NULL,
    species = "Orca",
    raw_blast_list = NULL
  ))

  testthat::expect_error(parsed_output <- parse_BLAST_json(
    filename = good_filename,
    species = 1,
    raw_blast_list = NULL
  ))

})
