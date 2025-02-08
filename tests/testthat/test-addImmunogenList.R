library(immunogenViewer)
exampleDF <- immunogenViewer::getProteinFeatures("P55087")
immunogenDF <- data.frame(
  start = c(10, 40, NA),
  end = c(30, 60, NA),
  seq = c(NA, NA, "RFKEAFSKAAQQTKGSYMEVEDNRSQVETDD"),
  name = c("A12", "B34", "HPA"),
  stringsAsFactors = FALSE
  )


test_that("Dataframe is returned if successul", {
  expect_s3_class(addImmunogenList(exampleDF, immunogenDF), "data.frame")
})

test_that("Missing column raises errors", {
  expect_error(addImmunogenList(exampleDF, data.frame(seq=c("RFKEAFSKAAQQTKGSYMEVEDNRSQVETDD"))),
    "The immunogen dataframe must contain a 'name' column.")
})

test_that("Missing column raises errors", {
  expect_error(addImmunogenList(exampleDF, data.frame(name=c("A12"))),
     "The immunogen dataframe must contain either a 'seq' column or both 'start' and 'end' columns.")
})

test_that("Wrong immunogen dataframe raises errors", {
  expect_error(addImmunogenList(exampleDF, data.frame(start = c(NA), end = c(NA), seq = c(NA), name = c("A12"))),
      "Invalid values. Please check your immunogen dataframe.")
})
