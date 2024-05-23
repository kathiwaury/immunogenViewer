test_that("Only strings are accepted as input", {
  expect_error(getProteinFeatures(NULL))
  expect_error(getProteinFeatures(4))
  expect_error(getProteinFeatures(3.14))
})

test_that("Dataframe is returned if successul", {
  expect_s3_class(getProteinFeatures("P55087"), "data.frame")
})


test_that("Only a list is accepted as input", {
  expect_error(getProteinSequence(2))
})


test_that("Protein sequence has to be a string", {
  expect_error(getProteinSequence(data.frame()))
})


test_that("Invalid URL throws error", {
  expect_error(accessUniprot("NotValid"))
  expect_error(accessPredictProtein("NotValid"))

})
