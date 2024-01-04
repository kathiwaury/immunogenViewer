#### add functionality: feature columns cannot be removed

#' Title
#'
#' @param proteinDF
#' @param name
#'
#' @return
#' @export
#'
#' @examples
removeImmunogen <- function(proteinDF, name) {

  colnamesDF <- colnames(proteinDF)

  if (checkIfImmunogenExists(colnamesDF, name)) {
    # proteinDF <- subset(proteinDF, select = -c(name))
    proteinDF <- proteinDF[,!(names(proteinDF) %in% name)]

    print(paste0("Successfully removed immunogen '", name, "'"))
    return(proteinDF)
  }

}


checkIfImmunogenExists <- function(colnamesDF, name) {

  if (name %in% colnamesDF) {
    return(TRUE)
  } else {
    stop("Immunogen not found in dataframe.")
  }

}
