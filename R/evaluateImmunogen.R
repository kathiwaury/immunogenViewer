# x <- getProteinFeatures("P55087")
# x <- addImmunogen(x, 0, 30, "A")
# x <- addImmunogen(x, 200, 250, "B")


#' Title
#'
#' @param proteinDF
#' @param immunogen
#'
#' @return
#' @export
#'
#' @examples
evaluateImmunogen <- function(proteinDF, immunogen) {

  features <- c("Uniprot", "Position", "Residue", "PTM", "disulfideBridge", "Membrane",
                "Binding", "Disorder", "secondaryStructure", "solventAccessibility")

  #### add functionality: feature columns cannot be used

  if (checkIfImmunogenExists(colnames(proteinDF), immunogen)) {

    immunogenDF <- filterProteinDataFrame(proteinDF, immunogen)

    print(paste0("Immunogen name: ", immunogen))

    start <- min(immunogenDF$Position)
    end <- max(immunogenDF$Position)
    seq <- paste(immunogenDF$Residue, collapse = "")

    print(paste0("Immunogen sequence: ", seq, " (Residues ", start, " - ", end, ")"))

    summaryDF <- createSummaryDataFrame(immunogenDF, immunogen)

    return(summaryDF)

  }
}


checkIfImmunogenExists <- function(colnamesDF, name) {

  if (name %in% colnamesDF) {
    return(TRUE)
  } else {
    stop("Immunogen not found in dataframe.")
  }

}


filterProteinDataFrame <- function(proteinDF, immunogen) {

  immunogenDF <- proteinDF[proteinDF[[immunogen]] == 1, ]

  return(immunogenDF)
}


calculateProportions <- function(immunogenDF, column) {

  proportions <- prop.table(table(immunogenDF[[column]]))

  return(proportions)
}


calculateRegionProportions <- function(immunogenDF, column) {

  if (sum(immunogenDF[[column]]) == 0) {
    return(0)

  } else {

    proportion <- prop.table(table(immunogenDF[[column]]))[["1"]]
    return(proportion)
  }

}


addMissingClasses <- function(proportionTable, classes) {

  for (i in classes) {

    if (!(i %in% names(proportionTable))) {
      proportionTable[[i]] = 0
    }
  }

  print(proportionTable)

  return(proportionTable)
}


createSummaryDataFrame <- function(immunogenDF, immunogen) {

  sumPTM <- sum(immunogenDF[["PTM"]])
  sumBridge <- sum(immunogenDF[["disulfideBridges"]])

  proportionMembrane <- calculateRegionProportions(immunogenDF, "Membrane")
  proportionDisorder <- calculateRegionProportions(immunogenDF, "Disorder")
  proportionBinding <- calculateRegionProportions(immunogenDF, "Binding")

  proportionsSecondaryStr <- prop.table(table(immunogenDF[["secondaryStructure"]]))
  proportionsSolventAcc <- prop.table(table(immunogenDF[["solventAccessibility"]]))

  # if values are not present in immunogen, add them as 0
  proportionsSecondaryStr <- addMissingClasses(proportionsSecondaryStr, c("Helix", "Sheet", "Other"))
  proportionsSolventAcc <- addMissingClasses(proportionsSolventAcc, c("Buried", "Exposed"))

  # create a summary dataframe
  summaryDF <- data.frame(
    Sum_PTM = sumPTM,
    Sum_DisulfideBridges = sumBridge,
    Proportion_Membrane = proportionMembrane,
    Proportion_Disorder = proportionDisorder,
    Proportion_Binding = proportionBinding,
    Proportions_Helix = proportionsSecondaryStr["Helix"][[1]],
    Proportions_Sheet = proportionsSecondaryStr["Sheet"][[1]],
    Proportions_Coil = proportionsSecondaryStr["Other"][[1]],
    Proportions_SolventAccessibility_Buried = proportionsSolventAcc["Buried"][[1]],
    Proportions_SolventAccessibility_Exposed = proportionsSolventAcc["Exposed"][[1]]
    )

  rownames(summaryDF) <- c(immunogen)

  return(summaryDF)
}
