
test_that("Checking valid user input for outputting predictions", {

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  #testthat::expect_equal(output_predictions(
  #  predictions = predictions
  #), NULL)
  # commented out for now, because it keeps printing the sequences out
  # (the test passes, though)
})

test_that("Checking for invalid user input for outputting predictions", {

  testthat::expect_error(output_predictions(
    predictions = NULL
  ))

  testthat::expect_error(output_predictions(
    predictions = 1
  ))

  testthat::expect_error(output_predictions(
    predictions = list()
  ))
})

# tests on file creation later on
