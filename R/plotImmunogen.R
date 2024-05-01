#' Plot protein features of one immunogen region
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param immunogen String, identifier name of immunogen
#'
#' @description
#' `plotImmunogen()` creates multiple ggplot objects within one figure. An Immunogen DataFrame is created by
#' filtering the Protein DataFrame for the relevant immunogen segment. For each feature with annotations in
#' the Immunogen DataFrame a plot is created. On the x axis the amino acid sequence of the immunogen is shown.
#'
#'
#' @return A ggplot object
#' @export
#'
#' @import ggplot2
#' @import patchwork
#'
#' @examples
#' proteinDF <- getProteinFeatures("P55087")
#' proteinDF <- addImmunogen(proteinDF, start=10, end=30, name="A12")
#' plotImmunogen(exampleDF, "A12")
plotImmunogen <- function(proteinDF, immunogen) {

  # check for valid immunogen names
  checkIfImmunogenExists(colnames(proteinDF), immunogen)

  # filter for relevant rows in protein dataframe
  immunogenDF <- filterProteinDataFrame(proteinDF, immunogen)

  # suppress warnings about missing values
  suppressWarnings({
    # plot immunogen features, overwrite plot limits to visualize full immunogen sequence
    structurePlot <- updatePlotLimits(plotSecondaryStructure(immunogenDF, NULL), immunogenDF)
    accessPlot <- updatePlotLimits(plotAccessibility(immunogenDF, NULL), immunogenDF)
    membranePlot <- updatePlotLimits(plotRegions(immunogenDF, NULL, "Membrane", "navyblue"), immunogenDF)
    bindingPlot <- updatePlotLimits(plotRegions(immunogenDF, NULL, "Binding", "forestgreen"), immunogenDF)
    disorderPlot <- updatePlotLimits(plotRegions(immunogenDF, NULL, "Disorder", "gold"), immunogenDF)
    ptmPlot <- updatePlotLimits(plotSinglePositions(immunogenDF, NULL, "PTM", "maroon"), immunogenDF)
    bridgesPlot <- updatePlotLimits(plotSinglePositions(immunogenDF, NULL, "disulfideBridge", "cornsilk4"), immunogenDF)

    # arrange plots below each other
    # print statement necessary to suppress warnings
    print(structurePlot / accessPlot  / membranePlot / bindingPlot / disorderPlot / ptmPlot / bridgesPlot)
  })
}


updatePlotLimits <- function(immunogenPlot, immunogenDF) {

  # suppress warning for already present x scale
  # expand_limits necessary to plot full immunogen sequence after setting labels to residues
  suppressMessages({
    immunogenPlot <- immunogenPlot +
    scale_x_continuous(breaks = min(immunogenDF$Position):max(immunogenDF$Position), labels = immunogenDF$Residue) +
    expand_limits(x = c(min(immunogenDF$Position), max(immunogenDF$Position)))
    })

  return(immunogenPlot)

}

