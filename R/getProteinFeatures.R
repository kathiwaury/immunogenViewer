# x <- getProteinFeatures("P55087")

#' Title
#'
#' @param uniprot
#'
#' @return
#' @export
#'
#' @examples
getProteinFeatures <- function(uniprot) {

  if (is.character(uniprot) == FALSE) {
    stop("Please provide a UniProt ID.")
  }

  resultUniprot <- accessUniprot(uniprot)
  uniprotDF <- createUniprotDataFrame(resultUniprot)

  predictProteinDF <- accessPredictProtein(uniprot)

  seq <- getProteinSequence(resultUniprot)
  df <- createFeatureDataFrame(uniprot, seq, uniprotDF, predictProteinDF)

  return(df)
}


accessUniprot <- function(uniprot) {

  # create protein-specific UniProt URL
  url <- paste0("https://www.uniprot.org/uniprot/", uniprot, ".json")

  # make GET request
  response <- httr::GET(url = url)

  # check if request was successful
  if (response$status_code == 200) {

    # parse JSON response
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding="UTF-8"))
    return(result)

  } else {
    stop(c("Error fetching data from Uniprot for ", uniprot))
  }
}


createUniprotDataFrame <- function(result) {

  if (is.list(result) == FALSE) {
    stop("Expected a list.")
  }

  uniprotDF <- result$features

  # any filtering to be added?

  return(uniprotDF)
}


accessPredictProtein <- function(uniprot) {

  # create protein-specific URL
  url <- paste0("https://api.predictprotein.org/v1/results/", uniprot)

  # make GET request
  response <- httr::GET(url = url)

  # check if the request was successful
  if (response$status_code == 200) {

    # parse JSON response
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding="UTF-8"))
    return(result$features)

  } else {
    stop(c("Error fetching data from PredictProtein for ", uniprot))
  }
}


getProteinSequence <- function(result) {

  if (is.list(result) == FALSE) {
    stop("Expected a list.")
  }

  seq <- as.character(result$sequence$value)

  if (is.character(seq) == FALSE || length(seq) == 0) {
    stop("Expected a character.")
  }

  return(seq)
}


createFeatureDataFrame <- function(uniprot, seq, uniprotDF, predictProteinDF) {

  # create vectors of positions and residues
  uniprotVector <- rep(uniprot, nchar(seq))
  numVector <- seq(1, nchar(seq))
  resVector <- strsplit(seq, NULL)[[1]]

  # create a dataframe with columns for positions and residues
  df <- data.frame(Uniprot = uniprotVector, Position = numVector, Residue = resVector)

  # PTMs (UniProt)
  PTMVector <- retrievePTMPositions(uniprotDF)
  df <- addPositions(df, PTMVector, "PTM")

  # disulfide bridges (UniProt)
  disulfideBridgeVector <- retrieveDisulfideBridgePositions(uniprotDF)
  df <- addPositions(df, disulfideBridgeVector, "disulfideBridge")

  # membrane residues  (UniProt)
  membraneVector <- retrieveMembranePositions(uniprotDF)
  df <- addPositions(df, membraneVector, "Membrane")

  # structure (PredictProtein)
  secStrVector <- createFeatureVector(predictProteinDF, "SECONDARY_STRUCTURE_(REPROF)")
  df["secondaryStructure"] <- secStrVector

  # solvent accessibility (PredictProtein)
  solAccVector <- createFeatureVector(predictProteinDF, "SOLVENT_ACCESSIBILITY_(REPROF)")
  df["solventAccessibility"] <- solAccVector

  # disorder (PredictProtein)
  disorderVector <- retrieveDisorderPositions(predictProteinDF)
  df <- addPositions(df, disorderVector, "Disorder")

  # binding (PredictProtein)
  bindingVector <- retrieveBindingPositions(predictProteinDF)
  df <- addPositions(df, bindingVector, "Binding")

  return(df)
}


addPositions <- function(df, featureVector, columnName) {

  df[[columnName]] <- ifelse(df$Position %in% featureVector, 1, 0)

  return(df)
}


collectRegionPositions <- function(start, end) {

  return(seq(start, end))
}


retrieveMembranePositions <- function(uniprotDF) {

  filters <- c("Transmembrane", "Intramembrane")
  filteredDF <- uniprotDF[uniprotDF$type %in% filters, ]

  membraneVector <- sort(unlist(mapply(collectRegionPositions, filteredDF$location$start$value,
      filteredDF$location$end$value, SIMPLIFY = FALSE)))

  return(membraneVector)
}


retrievePTMPositions <- function(uniprotDF) {

  filters <- c("Modified residue", "Lipidation", "Glycosylation")
  filteredDF <- uniprotDF[uniprotDF$type %in% filters, ]

  PTMVector <- sort(filteredDF$location$start$value)

  return(PTMVector)
}


retrieveDisulfideBridgePositions <- function(uniprotDF) {

  filters <- c("Disulfide bond")
  filteredDF <- uniprotDF[uniprotDF$type %in% filters, ]

  # disulfide bridges are reported as pairs, both begin and end position are a cysteine
  disulfideBridgeVector <- c(filteredDF$location$start$value, filteredDF$location$end$value)

  return(disulfideBridgeVector)
}


retrieveBindingPositions <- function(predictProteinDF) {

  # could add DNA and RNA binding as well
  filters <- c("PROTEIN_BINDING_(PRONA)")
  filteredDF <- predictProteinDF[predictProteinDF$type %in% filters, ]

  bindingVector <- sort(unlist(mapply(collectRegionPositions, filteredDF$begin,
      filteredDF$end, SIMPLIFY = FALSE)))

  return(bindingVector)
}


retrieveDisorderPositions <- function(predictProteinDF) {

  filters <- c("DISORDERED_REGION_(META-DISORDER)")
  filteredDF <- predictProteinDF[predictProteinDF$type %in% filters, ]

  disorderVector <- sort(unlist(mapply(collectRegionPositions, filteredDF$begin,
      filteredDF$end, SIMPLIFY = FALSE)))

  return(disorderVector)
}


createFeatureVector <- function(df, type) {

  filteredDF <- df[df$type == type, ]

  # use mapply with rep to generate the vector for each row
  # returns a named list that is not needed, unname() converts to simple vector
  featureVector <- unname(unlist(mapply(rep, filteredDF$description, filteredDF$end - filteredDF$begin + 1)))

  return(featureVector)
}
