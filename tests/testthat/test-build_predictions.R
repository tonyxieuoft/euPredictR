library(euPredictR)


test_that("Checking for existence of predictions for build_predictions provided
          valid input", {

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  testthat::expect_s3_class(predictions, "GenePredictions")
  testthat::expect_length(names(predictions), 5)
  testthat::expect_length(names(predictions$T_truncatus), 5)
  testthat::expect_type(predictions$O_orca$IL10, "character")
})


test_that("Checking for invalid user input for build_predictions", {

  testthat::expect_error(predictions <- build_predictions(
    raw_blast_list = NULL
  ))

  testthat::expect_error(predictions <- build_predictions(
    raw_blast_list = 1
  ))

  testthat::expect_error(predictions <- build_predictions(
    raw_blast_list = list()
  ))


})
