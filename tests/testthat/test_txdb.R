
test_that(desc = "txdb_to_gtf throws expected errors", {
  expect_error(txdb_to_gtf(txdb = NULL, file = "annotations.gtf"), 
               regexp = "class(txdb) not equal to \"TxDb\"", fixed = TRUE)
})

