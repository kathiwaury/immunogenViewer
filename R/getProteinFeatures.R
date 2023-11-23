getProteinFeatures <- function(uniprot) {
  print("Hello, world!")
  return(NULL)
}


accessUniprot <- function(uniprot) {

  if (is.character(uniprot) == FALSE) {
    warning("Please provide a valid UniProt ID.")
    return(NULL)
  }

  # create protein-specific UniProt URL
  url <- paste0("https://www.uniprot.org/uniprot/", uniprot, ".json")

  # make the GET request
  response <- httr::GET(url = url)

  print(response$status_code == 200)

  # check if request was successful
  if (response$status_code == 200) {

    # parse JSON response
    result <- jsonlite::fromJSON(httr::content(response, "text", encoding="UTF-8"))
    return(result)

  } else {

    warning("Error fetching data from UniProt.")
    return(NULL)
  }
}


# accessPredictProtein <- function(uniprot) {
#
#   # create protein-specific UniProt URL
#   url <- paste0("https://api.predictprotein.org/v1/results/", uniprot)
#
#   # Make the GET request
#   response <- http::GET(url = url)
#
#   # Check if the request was successful
#   if (http_status(response)$category == "Success") {
#
#     # Parse JSON response
#     result <- fromJSON(content(response, "text"))
#     return(result)
#
#   } else {
#
#     warning("Error fetching data from UniProt.")
#     return(NULL)
#   }
# }
