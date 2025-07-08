#' Plot protein features of one immunogen region
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param immunogen String, identifier name of immunogen
#'
#' @description
#' `plotImmunogen()` creates multiple ggplot objects within one figure. An Immunogen DataFrame is created by
#' filtering the Protein DataFrame for the relevant immunogen segment. A plot is created for each feature with annotations in
#' the Immunogen DataFrame. The amino acid sequence of the immunogen is shown on the x axis.
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
#' plotImmunogen(proteinDF, "A12")
plotImmunogen <- function(proteinDF, immunogen) {

  # check for valid immunogen names
  checkIfImmunogenExists(colnames(proteinDF), immunogen)

  # filter for relevant rows in protein dataframe
  immunogenDF <- filterProteinDataFrame(proteinDF, immunogen)

  # suppress warnings about missing values
  suppressWarnings({
    # plot immunogen features, overwrite plot limits to visualize full immunogen sequence
    structurePlot <- updatePlotLimits(plotSecondaryStructure(immunogenDF, NULL), immunogenDF)
    # accessPlot <- updatePlotLimits(plotAccessibility(immunogenDF, NULL), immunogenDF)
    membranePlot <- updatePlotLimits(plotRegions(immunogenDF, NULL, "Membrane", "navyblue",
        "Membrane (UniProt)"), immunogenDF)
    bindingPlot <- updatePlotLimits(plotRegions(immunogenDF, NULL, "ProteinBinding", "forestgreen",
        "Protein Binding (PredictProtein)"), immunogenDF)
    # disorderPlot <- updatePlotLimits(plotRegions(immunogenDF, NULL, "Disorder", "gold",
    #     "Disorder (PredictProtein)"), immunogenDF)
    ptmPlot <- updatePlotLimits(plotSinglePositions(immunogenDF, NULL, "PTM", "maroon",
        "PTM (UniProt)"), immunogenDF)
    bridgesPlot <- updatePlotLimits(plotSinglePositions(immunogenDF, NULL, "DisulfideBridge",
        "cornsilk4", "Disulfide Bridge (UniProt)"), immunogenDF)

    # arrange plots below each other
    # print statement necessary to suppress warnings
    # print(structurePlot / accessPlot  / membranePlot / bindingPlot / disorderPlot / ptmPlot / bridgesPlot)
    print(structurePlot / membranePlot / bindingPlot / ptmPlot / bridgesPlot)
  })
}


updatePlotLimits <- function(immunogenPlot, immunogenDF) {

  # suppress warning for already present x scale
  # expand_limits necessary to plot full immunogen sequence after setting labels to residues
  suppressMessages({
    immunogenPlot <- immunogenPlot +
    scale_x_continuous(breaks = min(immunogenDF$Position):max(immunogenDF$Position),
                       labels = immunogenDF$Residue) +
    expand_limits(x = c(min(immunogenDF$Position), max(immunogenDF$Position)))
    })

  return(immunogenPlot)

}

