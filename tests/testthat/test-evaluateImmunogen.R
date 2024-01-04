library(immunogenViewer)
exampleDF <- immunogenViewer::getProteinFeatures("P55087")
exampleDFWithImmunogen <- immunogenViewer::addImmunogen(exampleDF, 10, 30, "A12")


test_that("Dataframe is returned if successul", {
  expect_s3_class(evaluateImmunogen(exampleDFWithImmunogen, "A12"), "data.frame")
})

test_that("Non-exisiting immunogen name raises errors", {
  expect_error(evaluateImmunogen(exampleDF, "A12"), "Immunogen not found in dataframe.")
})
