library(immunogenViewer)
exampleDF_A <- immunogenViewer::getProteinFeatures("P55087") # no disulfide bridges (single position)
exampleDF_B <- immunogenViewer::getProteinFeatures("P07196") # no transmembrane (region)
exampleDFWithImmunogen <- immunogenViewer::addImmunogen(exampleDF_A, start=10, end=30, name="A12")
examplePlot_A <- immunogenViewer::plotImmunogen(exampleDFWithImmunogen, "A12")
# examplePlot_B <- immunogenViewer::plotProtein(exampleDF_B)
# examplePlotWithImmunogen <- immunogenViewer::plotProtein(exampleDFWithImmunogen)

test_that("ggplot is returned if successul", {
  expect_s3_class(examplePlot_A, "gg")
})

