#' Rename an existing immunogen
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param oldName String, current name of immunogen
#' @param newName String, new name of immunogen
#'
#' @return Updated Protein DataFrame with immunogen column renamed
#' @export
#'
#' @description
#' An existing immunogen is renamed in a Protein DataFrame by calling `renameImmunogen()`.
#'
#' @examples
#' proteinDF <- getProteinFeatures("P55087")
#' proteinDF <- addImmunogen(proteinDF, start=10, end=30, name="A12")
#' proteinDF <- renameImmunogen(proteinDF, "A12", "B12")
renameImmunogen <- function(proteinDF, oldName, newName) {

  # check for valid old immunogen names
  if (checkIfImmunogenExists(colnames(proteinDF), oldName)) {

    # check validity of new immunogen name
    colName <- checkImmunogenName(newName, colnames(proteinDF))

    # update column name
    names(proteinDF)[names(proteinDF) == oldName] <- newName

    return(proteinDF)
  }
}
