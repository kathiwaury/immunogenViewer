library(immunogenViewer)
exampleDF <- immunogenViewer::getProteinFeatures("P55087")
exampleDFWithImmunogen <- immunogenViewer::addImmunogen(exampleDF, start=200, end=230, name="A12")


test_that("Dataframe is returned if successul", {
  expect_s3_class(evaluateImmunogen(exampleDFWithImmunogen, "A12"), "data.frame")
  expect_s3_class(evaluateImmunogen(exampleDFWithImmunogen), "data.frame")
})

test_that("Non-exisiting immunogen name raises errors", {
  expect_error(evaluateImmunogen(exampleDF, "A12"), "Immunogen not found in dataframe.")
  expect_error(evaluateImmunogen(exampleDF, "PTM"), "The immunogen name is not valid.")
})
