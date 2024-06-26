library(immunogenViewer)
exampleDF <- immunogenViewer::getProteinFeatures("P55087")
colNames <- colnames(exampleDF)


test_that("Dataframe is returned if successul", {
  expect_s3_class(addImmunogen(exampleDF, start=10, end=30, name="A12"), "data.frame")
  expect_s3_class(addImmunogen(exampleDF, seq="RFKEAFSKAAQQTKGSYMEVEDNRSQVETDDLILKPGVVH", name="HPA"), "data.frame")
})

test_that("Invalid or missing sequence raises error", {
  expect_error(addImmunogen(exampleDF, seq="MFKEAFSKAAQQTKGSYMEVEDNRSQVEHDDLILKPGVPH", name="HPA"),
               "Immunogen sequence not found in the protein sequence.")
  expect_error(addImmunogen(exampleDF, name="A12"),
               "Either provide the immunogen sequence or start and end positions.")
})

test_that("Wrong immunogen range raises errors", {
  expect_error(checkImmunogenRange("A", "B", 100), "Immunogen start and end positions must be whole numbers.")
  expect_error(checkImmunogenRange(4.5, 9.333, 100), "Immunogen start and end positions must be whole numbers.")
  expect_error(checkImmunogenRange(10, 10, 100),
      "Immunogen start position must be smaller than immunogen end position.")
  expect_error(checkImmunogenRange(90, 110, 100), "Immunogen range must be within protein sequence length.")
  expect_error(checkImmunogenRange(20, 90, 100), "Immunogen should be between 10 and 50 residues long.")
})

test_that("Invalid immunogen name raises errors", {
  expect_error(checkImmunogenName("PTM", colNames),
      "Immunogen name cannot be same as existing features or immunogen column names.")
  expect_warning(checkImmunogenName(mean, colNames),
      "Unable to convert provided immunogen name to character. Using default name 'Immunogen' instead.")
  expect_warning(checkImmunogenName(NULL, colNames),
      "Unable to convert provided immunogen name to character. Using default name 'Immunogen' instead.")
})
