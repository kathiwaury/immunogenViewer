test_that("Only strings are accepted for UniProt IDs", {
  expect_null(accessUniprot(NULL))
  expect_null(accessUniprot(4))
  expect_null(accessUniprot(3.14))
})

test_that("Warning is raised if UniProt ID is incorrect", {
  expect_null(accessUniprot("Hello"))
  expect_null(accessUniprot("A00000"))
})
