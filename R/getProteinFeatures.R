#' Retrieve structural and functional features to create a protein DataFrame
#'
#' @param uniprot String, UniProt ID
#'
#' @description
#' By providing a valid UniProt ID, information from UniProt (https://www.uniprot.org/) and PredictProtein (https://predictprotein.org/) is queried via their
#' respective APIs. The retrieved information regarding secondary structure, solvent accessibility, membrane
#' regions, protein-binding regions, disordered regions, PTMs and disulfide bridges is saved per residue
#' within a Protein DataFrame. After calling `getProteinFeatures()`, immunogens can be added to the Protein
#' DataFrame.
#'
#' @return Protein DataFrame
#' @export
#'
#' @examples getProteinFeatures("P55087")
getProteinFeatures <- function(uniprot) {

  # check if Uniprot ID is valid
  if (is.character(uniprot) == FALSE) {
    stop("Please provide a UniProt ID.")
  }

  # retrieve results from UniProt
  resultUniprot <- accessUniprot(uniprot)
  uniprotDF <- resultUniprot$features

  # retrieve results from PredictProtein
  predictProteinDF <- accessPredictProtein(uniprot)

  # retrieve protein sequence
  seq <- getProteinSequence(resultUniprot)

  # create feature dataframe
  df <- createFeatureDataFrame(uniprot, seq, uniprotDF, predictProteinDF)
  print(paste0("Successfully created protein dataframe for UniProt ", uniprot, "."))

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
  proteinDF <- data.frame(Uniprot = uniprotVector, Position = numVector, Residue = resVector)

  # structure (PredictProtein)
  secStrVector <- createFeatureVector(predictProteinDF, "SECONDARY_STRUCTURE_(REPROF)")
  proteinDF["secondaryStructure"] <- secStrVector

  # solvent accessibility (PredictProtein)
  solAccVector <- createFeatureVector(predictProteinDF, "SOLVENT_ACCESSIBILITY_(REPROF)")
  proteinDF["solventAccessibility"] <- solAccVector

  # membrane residues  (UniProt)
  membraneVector <- retrieveMembranePositions(uniprotDF)
  proteinDF <- addPositions(proteinDF, membraneVector, "Membrane")

  # binding (PredictProtein)
  bindingVector <- retrieveRegionPositions(predictProteinDF, "PROTEIN_BINDING_(PRONA)") # add DNA and RNA binding as well?
  proteinDF <- addPositions(proteinDF, bindingVector, "Binding")

  # disorder (PredictProtein)
  disorderVector <- retrieveRegionPositions(predictProteinDF, "DISORDERED_REGION_(META-DISORDER)")
  proteinDF <- addPositions(proteinDF, disorderVector, "Disorder")

  # PTMs (UniProt)
  PTMVector <- retrieveSinglePositions(uniprotDF, c("Modified residue", "Lipidation", "Glycosylation"))
  proteinDF <- addPositions(proteinDF, PTMVector, "PTM")

  # disulfide bridges (UniProt)
  disulfideBridgeVector <- retrieveSinglePositions(uniprotDF, c("Disulfide bond"))
  proteinDF <- addPositions(proteinDF, disulfideBridgeVector, "disulfideBridge")

  return(proteinDF)
}


addPositions <- function(proteinDF, featureVector, columnName) {

  # for single positions with annotation add 1, otherwise 0
  proteinDF[[columnName]] <- ifelse(proteinDF$Position %in% featureVector, 1, 0)

  return(proteinDF)
}


collectRegionPositions <- function(start, end) {

  # get all positions within one region
  return(seq(start, end))
}


retrieveSinglePositions <- function(uniprotDF, columns) {

  # filter for relevant rows in Uniprot dataframe
  filteredDF <- uniprotDF[uniprotDF$type %in% columns, ]

  # collect all start and end positions of PTM and disulfide bridge annotations
  singlePositionsVector <- sort(c(filteredDF$location$start$value, filteredDF$location$end$value))

  return(singlePositionsVector)
}


retrieveMembranePositions <- function(uniprotDF) {

  # filter for relevant rows in Uniprot dataframe
  filters <- c("Transmembrane", "Intramembrane")
  filteredDF <- uniprotDF[uniprotDF$type %in% filters, ]

  # get all positions within membrane region
  membraneVector <- sort(unlist(mapply(collectRegionPositions, filteredDF$location$start$value,
      filteredDF$location$end$value, SIMPLIFY = FALSE)))

  return(membraneVector)
}


retrieveRegionPositions <- function(predictProteinDF, columns) {

  # filter for relevant rows in Uniprot dataframe
  filteredDF <- predictProteinDF[predictProteinDF$type %in% columns, ]

  # get all positions within region
  regionVector <- sort(unlist(mapply(collectRegionPositions, filteredDF$begin,
      filteredDF$end, SIMPLIFY = FALSE)))

  return(regionVector)
}


createFeatureVector <- function(predictProteinDF, type) {

  # filter for relevant rows in Uniprot dataframe
  filteredDF <- predictProteinDF[predictProteinDF$type == type, ]

  # use mapply with rep to generate the vector for each row
  # returns a named list that is not needed, unname() converts to simple vector
  featureVector <- unname(unlist(mapply(rep, filteredDF$description, filteredDF$end - filteredDF$begin + 1)))

  return(featureVector)
}
