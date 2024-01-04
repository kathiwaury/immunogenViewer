createReport <- function(proteinDF, immunogen) {

    # Open a PDF file
  pdf("output.pdf")

  # Your R plotting code or content generation here
  print(plotProtein(proteinDF))
  print(evaluateImmunogen(proteinDF, immunogen))

  # Close the PDF file
  dev.off()
  dev.off()
}



