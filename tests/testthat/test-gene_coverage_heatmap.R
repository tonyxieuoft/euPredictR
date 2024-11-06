

test_that("Checking that when valid user input is provided,
          gene_coverage_heatmap produces valid results", {

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )


  testthat::expect_failure(expect_equal(gene_coverage_heatmap(
    predictions = predictions
  ), NULL))

})


test_that("Checking for invalid user input for gene_coverage_heatmap", {

  testthat::expect_error(gene_coverage_heatmap(
    predictions = NULL
  ))

  testthat::expect_error(gene_coverage_heatmap(
    predictions = 1
  ))

  testthat::expect_error(gene_coverage_heatmap(
    predictions = list()
  ))



})
