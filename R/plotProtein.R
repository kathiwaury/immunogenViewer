# x <- getProteinFeatures("P55087")
# x <- addImmunogen(x, 5, 30, "A")
# x <- addImmunogen(x, 200, 250, "B")


#' Title
#'
#' @param proteinDF
#'
#' @return
#' @export
#'
#' @import ggplot2
#' @import patchwork
#'
#' @examples
plotProtein <- function(proteinDF) {

  immunogenRanges <- checkForImmunogens(proteinDF)

  # par(mfrow=c(7,1))

  ptmPlot <- plotSinglePositions(proteinDF, immunogenRanges, "PTM", "maroon")
  bridgesPlot <- plotSinglePositions(proteinDF, immunogenRanges, "disulfideBridge", "navyblue")
  membranePlot <- plotRegions(proteinDF, immunogenRanges, "Membrane", "cornsilk4")
  bindingPlot <- plotRegions(proteinDF, immunogenRanges, "Binding", "gold")
  disorderPlot <- plotRegions(proteinDF, immunogenRanges, "Disorder", "forestgreen")
  structurePlot <- plotSecondaryStructure(proteinDF, immunogenRanges)
  accessPlot <- plotAccessibility(proteinDF, immunogenRanges)

  # plotSinglePositions(proteinDF, immunogenRanges, "PTM", "maroon")
  # plotSinglePositions(proteinDF, immunogenRanges, "disulfideBridge", "navyblue")
  # plotRegions(proteinDF, immunogenRanges, "Membrane", "cornsilk4")
  # plotRegions(proteinDF, immunogenRanges, "Binding", "gold")
  # plotRegions(proteinDF, immunogenRanges, "Disorder", "forestgreen")
  # plotSecondaryStructure(proteinDF, immunogenRanges)
  # plotAccessibility(proteinDF, immunogenRanges)

  # gridExtra::grid.arrange(ptmPlot, bridgesPlot, membranePlot, bindingPlot, disorderPlot,
  #              structurePlot, accessPlot, ncol = 1)

  # par(mfrow = c(1, 1))

  ptmPlot / bridgesPlot / membranePlot / bindingPlot / disorderPlot / structurePlot / accessPlot
}


plotSinglePositions <- function(proteinDF, immunogenRanges, column, color) {

  if (sum(proteinDF[[column]]) == 0) {

    return(NULL)

  } else {

    locations <- data.frame(Location = proteinDF$Position[proteinDF[[column]] == 1], y = 0.5)

    plot <- ggplot() +
            geom_point(data = locations, aes(x = Location, y = y), colour = color, size = 3)

    plot <- setPlotAppearance(plot, column, nrow(proteinDF))

    if (!is.null(immunogenRanges)) {

      plot <- addImmunogensToPlot(plot, immunogenRanges)

    }
    return(plot)
  }
}


plotRegions <- function(proteinDF, immunogenRanges, column, color) {

  if (sum(proteinDF[[column]]) == 0) {

    return(NULL)

  } else {

    plot <- ggplot(data = proteinDF) +
            geom_rect(data = proteinDF[proteinDF[[column]] == 1,],
                      aes(xmin = Position, xmax = Position + 1, ymin = 0, ymax = 1),
                      fill = color)

    plot <- setPlotAppearance(plot, column, nrow(proteinDF))

    if (!is.null(immunogenRanges)) {

      plot <- addImmunogensToPlot(plot, immunogenRanges)

    }

    return(plot)
  }
}


plotSecondaryStructure <- function(proteinDF, immunogenRanges) {

  colors <- c("Helix" = "lightblue", "Strand" = "lightgreen", "Other" = "lightcoral")

  plot <- ggplot(data = proteinDF,
                 aes(fill = secondaryStructure)) +
          geom_rect(data = proteinDF[proteinDF$secondaryStructure %in% names(colors),],
                    aes(xmin = Position, xmax = Position + 1, ymin = 0, ymax = 1)) +
          scale_fill_manual(values = colors)

  plot <- setPlotAppearance(plot, "Secondary Structure", nrow(proteinDF))

  if (!is.null(immunogenRanges)) {

    plot <- addImmunogensToPlot(plot, immunogenRanges)

  }

  return(plot)

}


plotAccessibility <- function(proteinDF, immunogenRanges) {

  colors <- c("Buried" = "grey", "Exposed" = "aquamarine4")

  plot <- ggplot(data = proteinDF, ggplot2::aes(fill = solventAccessibility)) +
          geom_rect(data = proteinDF[proteinDF$solventAccessibility %in% names(colors),],
                    aes(xmin = Position, xmax = Position + 1, ymin = 0, ymax = 1)) +
          scale_fill_manual(values = colors)

  plot <- setPlotAppearance(plot, "Solvent Accessibility", nrow(proteinDF))

  if (!is.null(immunogenRanges)) {

    plot <- addImmunogensToPlot(plot, immunogenRanges)

  }

  return(plot)

}


checkForImmunogens <- function(proteinDF) {

  features <- c("Uniprot", "Position", "Residue", "PTM", "disulfideBridge", "Membrane",
                "Binding", "Disorder", "secondaryStructure", "solventAccessibility")

  if (length(setdiff(colnames(proteinDF), features)) == 0) {

    return(NULL)

  } else {

    immunogens <- names(proteinDF)[!(names(proteinDF) %in% features)]

    immunogenRanges <- list()

    for (i in immunogens) {

      start <- which(diff(c(0, proteinDF[[i]])) == 1)
      end <- which(diff(c(proteinDF[[i]], 0)) == -1)

      immunogenRanges[[i]] <- c(start, end)

    }
    return(immunogenRanges)
  }
}


addImmunogensToPlot <- function(plot, immunogenRanges) {
  for (i in names(immunogenRanges)) {
    rect_data <- data.frame(xmin = immunogenRanges[[i]][1],
                            xmax = immunogenRanges[[i]][2],
                            ymin = 0,
                            ymax = 1)

    plot <- plot +
      ggplot2::geom_rect(data = rect_data,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         fill = NA, color = "darkred", linewidth = 1.2) +
      ggplot2::geom_text(data = data.frame(x = mean(c(rect_data$xmin, rect_data$xmax)), y = 1, label = i),
                         aes(x = x, y = y, label = label), vjust = -0.5, hjust = 0.5, size = 3,
                         inherit.aes = FALSE)
  }

  return(plot)
}


setPlotAppearance <- function(plot, column, seqLength) {

  plot <- plot +
    theme_classic() +
    lims(x = c(0, seqLength), y = c(0, 2)) +
    labs(x = "Residue", title = column) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(size=10),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size=8))

  return(plot)
}
