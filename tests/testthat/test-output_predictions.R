library(euPredictR)

test_that("Checking for valid data frame format produced when user input is
          correct", {

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  output_df <- output_predictions_df(
    predictions = predictions
  )

  testthat::expect_true(is.character(output_df$Gene))
  testthat::expect_true(is.character(output_df$Species))
  testthat::expect_true(is.character(output_df$Sequence))

  testthat::expect_contains(output_df$Species, names(predictions))
  testthat::expect_contains(output_df$Gene, names(predictions$O_orca))

})

test_that("Checks that output_predictions_fasta returns NULL when user
          input is correct", {

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  base::sink(file = nullfile())
  output_fasta <- output_predictions_fasta(
    predictions = predictions
  )
  base::sink()

  testthat::expect_null(output_fasta)
})

test_that("Checking for invalid user input for outputting predictions in
          fasta", {

  testthat::expect_error(output_predictions_fasta(
    predictions = NULL
  ))

  testthat::expect_error(output_predictions_fasta(
    predictions = 1
  ))

  testthat::expect_error(output_predictions_fasta(
    predictions = list()
  ))

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  testthat::expect_error(output_predictions_fasta(
    predictions = raw_blast_list,
    output_file = 25
  ))
})

test_that("Checking for invalid user input for outputting predictions in df", {

  testthat::expect_error(output_predictions_df(
    predictions = NULL
  ))

  testthat::expect_error(output_predictions_df(
    predictions = 1
  ))

  testthat::expect_error(output_predictions_df(
    predictions = list()
  ))
})
