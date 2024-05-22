#' Create a summary DataFrame of the structural and functional properties of an immunogen
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param immunogen String, identifier name of immunogen (if not defined, all immunogens are evaluated)
#'
#' @return Summary DataFrame providing statistics on immunogen
#' @export
#'
#' @description
#' By calling `evaluateImmunogen()`, the immunogens associated with a Protein DataFrame can be evaluated regarding
#' their suitability fir antibody binding in natively folded proteins. By calling the function without specifying
#' an immunogen, all immunogens of the current protein will be evaluated. The summary DataFrame contains one row per
#' evaluated immunogen.
#'
#' @examples
#' proteinDF <- getProteinFeatures("P55087")
#' proteinDF <- addImmunogen(proteinDF, start=10, end=30, "A12")
#' proteinDF <- addImmunogen(proteinDF, seq="RFKEAFSKAAQQTKGSYMEVEDNRSQVETDD", name="HPA")
#' evaluateImmunogen(proteinDF, "A12")
#' evaluateImmunogen(proteinDF)
evaluateImmunogen <- function(proteinDF, immunogen=NULL) {

  # feature column names of protein dataframe
  features <- c("Uniprot", "Position", "Residue", "PTM", "disulfideBridge", "Membrane",
                "Binding", "Disorder", "secondaryStructure", "solventAccessibility")

  fullDF <- data.frame()

  # if no immunogen is specified, all immunogens are included
  if (is.null(immunogen)) {

    print("No immunogen specified, evaluating all immunogens.")

    # select all column names that are not features
    immunogens <- colnames(proteinDF)[!colnames(proteinDF) %in% features]

  } else {

    # if immunogen specified, create list of only this immunogen
    immunogens <- c(immunogen)
  }

  # loop through list of immunogens
  for (i in immunogens) {

    # check for valid immunogen names
    if (checkIfImmunogenExists(colnames(proteinDF), i)) {

      # filter for relevant rows in protein dataframe
      immunogenDF <- filterProteinDataFrame(proteinDF, i)

      print(paste0("Immunogen name: ", i))

      # retrieve and print immungen sequence and position
      start <- min(immunogenDF$Position)
      end <- max(immunogenDF$Position)
      seq <- paste(immunogenDF$Residue, collapse = "")
      print(paste0("Immunogen sequence: ", seq, " (Residues ", start, " - ", end, ")"))

      # summarize features, add each immunogen summary to final dataframe
      summaryDF <- createSummaryDataFrame(immunogenDF, i)
      fullDF <- rbind(fullDF, summaryDF)

    }
  }

  return(fullDF)
}


checkIfImmunogenExists <- function(colnamesDF, name) {

  # feature column names of protein dataframe
  features <- c("Uniprot", "Position", "Residue", "PTM", "disulfideBridge", "Membrane",
                "Binding", "Disorder", "secondaryStructure", "solventAccessibility")

  # raise error if immunogen names not present in dataframe or same as feature column name
  if (!(name %in% colnamesDF)) {
    stop("Immunogen not found in dataframe.")

  } else {
    if (name %in% features) {
      stop("The immunogen name is not valid.")

    } else {
      return(TRUE)
    }
  }
}


filterProteinDataFrame <- function(proteinDF, immunogen) {

  # filter for relevant rows of protein dataframe
  immunogenDF <- proteinDF[proteinDF[[immunogen]] == 1, ]

  return(immunogenDF)
}


calculateRegionProportions <- function(immunogenDF, column) {

  # if no annotations, return 0
  if (sum(immunogenDF[[column]]) == 0) {
    return(0)

  } else {

    # calculate proportion of positive annotations ("1") within immunogen
    proportion <- prop.table(table(immunogenDF[[column]]))[["1"]]
    return(proportion)
  }

}


addMissingClasses <- function(proportionTable, classes) {

  for (i in classes) {

    # if value is not present, set proportion to 0
    if (!(i %in% names(proportionTable))) {
      proportionTable[[i]] = 0
    }
  }

  return(proportionTable)
}


createSummaryDataFrame <- function(immunogenDF, immunogen) {

  # count number of positions within immunogen
  sumPTM <- sum(immunogenDF[["PTM"]])
  sumBridge <- sum(immunogenDF[["disulfideBridges"]])

  proportionMembrane <- calculateRegionProportions(immunogenDF, "Membrane")
  proportionDisorder <- calculateRegionProportions(immunogenDF, "Disorder")
  proportionBinding <- calculateRegionProportions(immunogenDF, "Binding")

  # add secondary structure proportions
  proportionsSecondaryStr <- prop.table(table(immunogenDF[["secondaryStructure"]]))
  proportionsSecondaryStr <- addMissingClasses(proportionsSecondaryStr, c("Helix", "Sheet", "Other"))

  # add buried/exposed proportions
  proportionsSolventAcc <- prop.table(table(immunogenDF[["solventAccessibility"]]))
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

  # set all rownames to immunogen name
  rownames(summaryDF) <- c(immunogen)

  return(summaryDF)
}
