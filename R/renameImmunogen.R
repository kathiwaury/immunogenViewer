#' Title
#'
#' @param proteinDF
#' @param oldName
#' @param newName
#'
#' @return
#' @export
#'
#' @examples
renameImmunogen <- function(proteinDF, oldName, newName) {

  colnamesDF <- colnames(proteinDF)

  if (checkIfImmunogenExists(colnamesDF, oldName)) {
    names(proteinDF)[names(proteinDF) == oldName] <- newName
    print(paste0("Successfully updated immunogen name '", oldName, "' to '", newName, "'"))
    return(proteinDF)
  }
}


checkIfImmunogenExists <- function(colnamesDF, oldName) {

  if (oldName %in% colnamesDF) {
    return(TRUE)
  } else {
    stop("Immunogen not found in dataframe.")
  }

}
