#' Add multiple immunogens to the Protein DataFrame
#'
#' @param proteinDF Protein DataFrame created by call to getProteinFeatures()
#' @param immunogenDF DataFrame where each row represents an immunogen.
#'        Must contain columns: `start` (integer) and `end` (integer) or `seq` (string), and `name` (string).
#'
#' @description
#' Calls `addImmunogen()` for each row in `immunogenDF` to add multiple immunogens
#' to the given `proteinDF`.
#'
#' @return Updated Protein DataFrame with all immunogens added as new columns
#' @export
#'
#' @examples
#' proteinDF <- getProteinFeatures("P55087")
#' immunogenDF <- data.frame(
#'   start = c(10, 40, NA),
#'   end = c(30, 60, NA),
#'   seq = c(NA, NA, "RFKEAFSKAAQQTKGSYMEVEDNRSQVETDD"),
#'   name = c("A12", "B34", "HPA"),
#'   stringsAsFactors = FALSE
#' )
#' proteinDF <- addImmunogenList(proteinDF, immunogenDF)
addImmunogenList <- function(proteinDF, immunogenDF) {

  # Check if required columns exist
  if (!"name" %in% colnames(immunogenDF)) {
    stop("The immunogen dataframe must contain a 'name' column.")
  }
  if (!("seq" %in% colnames(immunogenDF)) && !all(c("start", "end") %in% colnames(immunogenDF))) {
    stop("The immunogen dataframe must contain either a 'seq' column or both 'start' and 'end' columns.")
  }

  # Iterate through each row and add immunogen
  for (i in seq_len(nrow(immunogenDF))) {

    row <- immunogenDF[i, ]

    # Call to addImmunogen based on available data
    if (!is.na(row$start) && !is.na(row$end)) {
      proteinDF <- addImmunogen(proteinDF, start = row$start, end = row$end, name = row$name)
    } else if (!is.na(row$seq)) {
      proteinDF <- addImmunogen(proteinDF, seq = as.character(row$seq), name = row$name)
    } else {
      stop("Invalid values. Please check your immunogen dataframe.")
    }
  }

  return(proteinDF)
}
