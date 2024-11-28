library(euPredictR)

test_that("Checking that when valid user input is provided,
          display_phylogeny produces valid results", {

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  png(tempfile(fileext = ".png"))
  phylogeny <- display_phylogeny(predictions = predictions,
                                 gene_name = "RHO")
  dev.off()

  testthat::expect_true(!is.null(phylogeny))
})

test_that("Checking for invalid user input for display_phylogeny", {

  testthat::expect_error(display_phylogeny(
    predictions = NULL,
    gene_name = "RHO"
  ))

  testthat::expect_error(display_phylogeny(
    predictions = 1,
    gene_name = "RHO"
  ))

  testthat::expect_error(display_phylogeny(
    predictions = list(),
    gene_name = "RHO"
  ))

  predictions <- build_predictions(
    raw_blast_list = sample_raw_blast_list
  )

  testthat::expect_error(display_phylogeny(
    predictions = predictions,
    gene_name = 2
  ))

  testthat::expect_error(display_phylogeny(
    predictions = predictions,
    gene_name = "invalid_gene_name"
  ))

})
