#' Title
#'
#' @param proteinDF
#' @param start
#' @param end
#' @param name
#'
#' @return
#' @export
#'
#' @examples
addImmunogen <- function(proteinDF, start, end, name) {

  proteinLength <- nrow(proteinDF)
  colnamesDF <- colnames(proteinDF)

  # check validity of immunogen name and range input
  colName <- checkImmunogenName(name, colnamesDF)
  if (checkImmunogenRange(start, end, proteinLength)) {

    # add immunogen as column to protein dataframe
    immunogenVector <- seq(start, end)
    proteinDF[[colName]] <- ifelse(proteinDF$Position %in% immunogenVector, 1, 0)

    return(proteinDF)
  }

}


checkImmunogenName <- function(name, colNamesDF) {

  colName <- tryCatch(
    expr = as.character(name),
    error = function(err) {
      warning("Unable to convert provided immunogen name to character. Using default name 'Immunogen' instead.")
      return("Immunogen")
    }
  )

  if (length(colName) == 0) {
    warning("Unable to convert provided immunogen name to character. Using default name 'Immunogen' instead.")
    colName <- "Immunogen"
  }


  if (colName %in% colNamesDF) {
    stop("Immunogen name cannot be same as existing feature or immunogen column names.")
  }

  return(colName)
}


checkImmunogenRange <- function(start, end, length) {

  if (is.numeric(start) == FALSE || is.numeric(end) == FALSE) {
    stop("Immunogen start and end positions must be whole numbers.")

  } else if (start %% 1 != 0 || end %% 1 != 0) {
      stop("Immunogen start and end positions must be whole numbers.")

  } else if (start >= end) {
    stop("Immunogen start position must be smaller than immunogen end position.")

  } else if (start >= length ||  end > length) {
    stop("Immunogen range must be within protein sequence length.")
  }

  immunogenLength <- end - start

  if (immunogenLength < 10 || immunogenLength > 50) {
    stop(c("Immunogen should be between 10 and 50 residues long."))
  }

  return(TRUE)
}
