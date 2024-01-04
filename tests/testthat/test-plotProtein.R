library(immunogenViewer)
exampleDF_A <- immunogenViewer::getProteinFeatures("P55087") # no disulfide bridges (single position)
exampleDF_B <- immunogenViewer::getProteinFeatures("P07196") # no transmembrane (region)
exampleDFWithImmunogen <- immunogenViewer::addImmunogen(exampleDF_A, 10, 30, "A12")
examplePlot_A <- immunogenViewer::plotProtein(exampleDF_A)
examplePlot_B <- immunogenViewer::plotProtein(exampleDF_B)
examplePlotWithImmunogen <- immunogenViewer::plotProtein(exampleDFWithImmunogen)

test_that("ggplot is returned if successul", {
  expect_s3_class(examplePlot_A, "gg")
  expect_s3_class(examplePlot_B, "gg")
  expect_s3_class(examplePlotWithImmunogen, "gg")
})
