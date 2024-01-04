library(immunogenViewer)
exampleDF <- immunogenViewer::getProteinFeatures("P55087")
exampleDF <- immunogenViewer::addImmunogen(exampleDF, 10, 30, "A12")
colNames <- colnames(exampleDF)

test_that("Dataframe is returned if successul", {
  expect_s3_class(renameImmunogen(exampleDF, "A12", "B12"), "data.frame")
})

test_that("Exisiting immunogen name returns TRUE", {
  expect_true(checkIfImmunogenExists(colNames, "A12"))
})

test_that("Non-exisiting immunogen name raises errors", {
  expect_error(checkIfImmunogenExists(colNames, "Hello"), "Immunogen not found in dataframe.")
})
