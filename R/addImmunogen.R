#' Add an immunogen to the Protein DataFrame
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param start Integer, start position of immunogen
#' @param end Integer, end position of immunogen
#' @param seq String, immunogen amino acid sequence (must be substring of protein sequence)
#' @param name String, chosen name to identify immunogen
#'
#' @description
#' An immunogen is associated with a protein by adding a column to a Protein DataFrame.
#' The immunogen is specified by a unique name. Its position is defined by either the start
#' and end position within the protein sequence or by supplying the immunogen sequence which
#' must be a substring of the protein's sequence.
#'
#' @return Updated Protein DataFrame with immunogen added as a new column
#' @export
#'
#' @examples
#' proteinDF <- getProteinFeatures("P55087")
#' proteinDF <- addImmunogen(proteinDF, start=10, end=30, name="A12")
#' proteinDF <- addImmunogen(proteinDF, seq="RFKEAFSKAAQQTKGSYMEVEDNRSQVETDD", name="HPA")
addImmunogen <- function(proteinDF, start=NULL, end=NULL, seq=NULL, name) {

  # save protein length and sequence
  proteinLength <- nrow(proteinDF)
  proteinSeq <- paste(proteinDF$Residue, collapse = "")

  colnamesDF <- colnames(proteinDF)

  # check validity of immunogen name
  colName <- checkImmunogenName(name, colnamesDF)

  # if immunogen sequence is provided, find position within protein sequence
  if (!is.null(seq)) {

    # find immunogen sequence start position
    start <- regexpr(seq, proteinSeq)[[1]]

    if (start == -1) {
      stop("Immunogen sequence not found in the protein sequence.")
    }

    # calculate immunogen end position
    end <- start + nchar(seq) - 1

  # if start and end position are provided, check validity of range
  } else if (!is.null(start) & !is.null(end)) {
    checkImmunogenRange(start, end, proteinLength)

  } else {
    stop("Either provide the immunogen sequence or start and end positions.")
  }

  # add immunogen as column to protein dataframe
  immunogenVector <- seq(start, end)
  proteinDF[[colName]] <- ifelse(proteinDF$Position %in% immunogenVector, 1, 0)

  print(paste0("Successfully added immunogen '", colName, "' to the protein dataframe."))
  return(proteinDF)

}


checkImmunogenName <- function(name, colNamesDF) {

  # if immunogen name cannot be converted to string, use default name
  name <- tryCatch(
    expr = as.character(name),
    error = function(err) {
      warning("Unable to convert provided immunogen name to character. Using default name 'Immunogen' instead.")
      return("Immunogen")
    }
  )

  if (length(name) == 0) {
    warning("Unable to convert provided immunogen name to character. Using default name 'Immunogen' instead.")
    colName <- "Immunogen"
  }

  # compare to existing column names
  if (name %in% colNamesDF) {
    stop("Immunogen name cannot be same as existing features or immunogen column names.")
  }

  return(name)
}


checkImmunogenRange <- function(start, end, length) {

  # check that start and end are integers
  if (is.numeric(start) == FALSE || is.numeric(end) == FALSE) {
    stop("Immunogen start and end positions must be whole numbers.")

  } else if (start %% 1 != 0 || end %% 1 != 0) {
      stop("Immunogen start and end positions must be whole numbers.")

  # check for valid start and end positions
  } else if (start >= end) {
    stop("Immunogen start position must be smaller than immunogen end position.")

  } else if (start >= length ||  end > length) {
    stop("Immunogen range must be within protein sequence length.")
  }

  # check for allowed immunogen length
  immunogenLength <- end - start

  if (immunogenLength < 10 || immunogenLength > 50) {
    stop(c("Immunogen should be between 10 and 50 residues long."))
  }

  return(TRUE)
}
