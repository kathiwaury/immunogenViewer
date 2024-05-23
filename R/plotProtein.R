#' Plot protein features with immunogens highlighted
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#'
#' @description
#' A call to `plotProtein()` visualizes all relevant protein features within one figure along
#' the entire protein sequence. All immunogens associated with the protein are highlighted at their
#' position along the protein sequence by darkred boxes.
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
#' plotProtein(proteinDF)
plotProtein <- function(proteinDF) {

  # collect immunogen ranges to add boxes to plots
  immunogenRanges <- checkForImmunogens(proteinDF)

  # suppress warnings about missing values
  suppressWarnings({
    # plot immunogen features, overwrite plot limits to visualize full immunogen sequence
    structurePlot <- plotSecondaryStructure(proteinDF, immunogenRanges)
    accessPlot <- plotAccessibility(proteinDF, immunogenRanges)
    membranePlot <- plotRegions(proteinDF, immunogenRanges, "Membrane", "navyblue")
    bindingPlot <- plotRegions(proteinDF, immunogenRanges, "Binding", "forestgreen")
    disorderPlot <- plotRegions(proteinDF, immunogenRanges, "Disorder", "gold")
    ptmPlot <- plotSinglePositions(proteinDF, immunogenRanges, "PTM", "maroon")
    bridgesPlot <- plotSinglePositions(proteinDF, immunogenRanges, "disulfideBridge", "cornsilk4")

    # arrange plots below each other
    # print statement necessary to suppress warnings
    print(structurePlot / accessPlot / membranePlot / bindingPlot / disorderPlot / ptmPlot / bridgesPlot )
  })
}


checkForImmunogens <- function(proteinDF) {

  # feature column names of protein dataframe
  features <- c("Uniprot", "Position", "Residue", "PTM", "disulfideBridge", "Membrane",
                "Binding", "Disorder", "secondaryStructure", "solventAccessibility")

  # if no immunogens are present, return NULL
  if (length(setdiff(colnames(proteinDF), features)) == 0) {

    return(NULL)

  } else {

    # all non-feature columns are immunogen columns
    immunogens <- names(proteinDF)[!(names(proteinDF) %in% features)]

    immunogenRanges <- list()

    # for every immunogen add start and end position to list
    for (i in immunogens) {

      start <- which(diff(c(0, proteinDF[[i]])) == 1)
      end <- which(diff(c(proteinDF[[i]], 0)) == -1)

      immunogenRanges[[i]] <- c(start, end)

    }

    return(immunogenRanges)
  }
}


plotSecondaryStructure <- function(proteinDF, immunogenRanges) {

  # set color for secondary structure values
  colors <- c("Helix" = "lightblue", "Strand" = "lightgreen", "Other" = "lightcoral")

  # plot rectangles for every block of secondary structure
  plot <- ggplot(data = proteinDF, aes(fill = secondaryStructure), alpha = 0.05) +
    geom_rect(data = proteinDF[proteinDF$secondaryStructure %in% names(colors),],
              aes(xmin = Position, xmax = Position + 1, ymin = 0, ymax = 1,
                  fill = secondaryStructure, color = secondaryStructure)) +
    # drop title of legend
    scale_fill_manual(values = colors, name = NULL) +
    # suppress legend to not display twice
    scale_color_manual(values = colors, guide = FALSE)

  plot <- setPlotAppearance(plot, "Secondary Structure", nrow(proteinDF))

  # if immunogens are present in protein dataframe, visualize in plot
  if (!is.null(immunogenRanges)) {

    plot <- addImmunogensToPlot(plot, immunogenRanges)

  }

  return(plot)

}


plotAccessibility <- function(proteinDF, immunogenRanges) {

  # set color for buried/exposed values
  colors <- c("Buried" = "grey", "Exposed" = "aquamarine4")

  # plot rectangles for every block of surface accessibility
  plot <- ggplot(data = proteinDF, aes(fill = solventAccessibility)) +
    geom_rect(data = proteinDF[proteinDF$solventAccessibility %in% names(colors),],
              aes(xmin = Position, xmax = Position + 1, ymin = 0, ymax = 1,
                  fill = solventAccessibility, color = solventAccessibility)) +
    # drop title of legend
    scale_fill_manual(values = colors, name = NULL) +
    # suppress legend to not display twice
    scale_color_manual(values = colors, guide = FALSE)

  plot <- setPlotAppearance(plot, "Solvent Accessibility", nrow(proteinDF))

  # if immunogens are present in protein dataframe, visualize in plot
  if (!is.null(immunogenRanges)) {

    plot <- addImmunogensToPlot(plot, immunogenRanges)

  }

  return(plot)

}


plotRegions <- function(proteinDF, immunogenRanges, column, color) {

  # if no annotations, do not plot anything
  if (sum(proteinDF[[column]]) == 0) {

    return(NULL)

  } else {

    # plot rectangles for every region
    plot <- ggplot(data = proteinDF) +
      geom_rect(data = proteinDF[proteinDF[[column]] == 1,],
                aes(xmin = Position, xmax = Position + 1, ymin = 0, ymax = 1),
                fill = color, color = color)

    plot <- setPlotAppearance(plot, column, nrow(proteinDF))

    # if immunogens are present in protein dataframe, visualize in plot
    if (!is.null(immunogenRanges)) {

      plot <- addImmunogensToPlot(plot, immunogenRanges)

    }

    return(plot)
  }
}


plotSinglePositions <- function(proteinDF, immunogenRanges, column, color) {

  # if no annotations, do not plot anything
  if (sum(proteinDF[[column]]) == 0) {

    return(NULL)

  } else {

    locations <- data.frame(Location = proteinDF$Position[proteinDF[[column]] == 1], y = 0.5)

    # plot points for every position
    plot <- ggplot() +
            geom_point(data = locations, aes(x = Location, y = y), colour = color, size = 3)

    plot <- setPlotAppearance(plot, column, nrow(proteinDF))

    # if immunogens are present in protein dataframe, visualize in plot
    if (!is.null(immunogenRanges)) {

      plot <- addImmunogensToPlot(plot, immunogenRanges)

    }

    return(plot)
  }
}


addImmunogensToPlot <- function(plot, immunogenRanges) {

  # create borders for every immunogen range
  for (i in names(immunogenRanges)) {
    rect_data <- data.frame(xmin = immunogenRanges[[i]][1],
                            xmax = immunogenRanges[[i]][2],
                            ymin = 0,
                            ymax = 1)

    # plot dark red rectangles around feature plots for every immunogen
    plot <- plot +
      ggplot2::geom_rect(data = rect_data,
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                         fill = NA, color = "darkred", linewidth = 1.2) +
      # add name of immunogen above box
      ggplot2::geom_text(data = data.frame(x = mean(c(rect_data$xmin, rect_data$xmax)), y = 1, label = i),
                         aes(x = x, y = y, label = label), vjust = -0.5, hjust = 0.5, size = 3,
                         inherit.aes = FALSE)
  }

  return(plot)
}


setPlotAppearance <- function(plot, column, seqLength) {

  # update plot settings
  plot <- plot +
    theme_classic() +
    lims(x = c(0, seqLength), y = c(0, 2)) +
    labs(x = "Residue", title = column) +
    # do not show y axis
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          plot.title = element_text(size=10),
          axis.title.x = element_text(size=8),
          axis.text.x = element_text(size=8))

  return(plot)
}
