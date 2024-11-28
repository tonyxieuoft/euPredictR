#context("Test the parse_BLAST_json file in ParseInput.R")
library(euPredictR)


test_that("Valid .json file for dolphin and orca IL6,IL10,TNF,RHO, SOX2 HSPs
          obtained by BLAST (see data.R for details)", {
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
     species = "O_orca",
     raw_blast_list = raw_blast_list1
   )

   testthat::expect_s3_class(raw_blast_list1, "RawBlastList")
   testthat::expect_s3_class(raw_blast_list2, "RawBlastList")
})


test_that("Valid directory containing .json files for dolphin, orca, mouse,
          lamprey, and zebrafish IL6,IL10,TNF,RHO, SOX2 HSPs obtained by BLAST
          (see data.R for details)", {
            good_dirname <- system.file("extdata/example_raw_output_directory",
                                          package = "euPredictR")

            raw_blast_list <- parse_multiple_BLAST_json(
              dir_path = good_dirname,
              raw_blast_list = NULL
            )

            species_names <- names(raw_blast_list)

            testthat::expect_s3_class(raw_blast_list, "RawBlastList")
            testthat::expect_equal(length(species_names), 5)

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

test_that("Check for valid .json but NOT BLAST format file for
          parse_BLAST_json. Specifically, the file misses a single
          query_title entry.", {

  good_filename <- system.file("extdata",
                               "invalid_file_non_BLAST_format.json",
                               package = "euPredictR")
  testthat::expect_error(parsed_output <- parse_BLAST_json(
    filename = good_filename,
    species = "Orca",
    raw_blast_list = NULL
  ))
})


test_that("Checking for invalid user input for parse_multiple_BLAST_json. More
          specifically, empty directory inputs and input directories that
          contain non .json files", {

  bad_dir1 <- system.file("extdata/invalid_directory_for_test1",
                               package = "euPredictR")
  bad_dir2 <- system.file("extdata/invalid_directory_for_test2",
                          package = "euPredictR")

  testthat::expect_error(parsed_output <- parse_multiple_BLAST_json(
    dir_path = bad_dir1,
    raw_blast_list = NULL
  ))

  testthat::expect_error(parsed_output <- parse_multiple_BLAST_json(
    dir_path = bad_dir2,
    raw_blast_list = NULL
  ))

})
