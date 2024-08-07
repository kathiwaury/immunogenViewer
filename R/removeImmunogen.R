#' Remove an existing immunogen
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param name String, name of immunogen
#'
#' @return Updated Protein DataFrame with immunogen column removed
#' @export
#'
#' @description
#' An existing immunogen is removed from a Protein DataFrame by calling `removeImmunogen()`.
#'
#'
#' @examples
#' proteinDF <- getProteinFeatures("P55087")
#' proteinDF <- addImmunogen(proteinDF, start=10, end=30, name="A12")
#' proteinDF <- removeImmunogen(proteinDF, "A12")
removeImmunogen <- function(proteinDF, name) {

  # check for valid immunogen names
  if (checkIfImmunogenExists(colnames(proteinDF), name)) {

    # remove immunogen column from protein dataframe
    proteinDF <- proteinDF[,!(names(proteinDF) %in% name)]

    return(proteinDF)
  }

}
