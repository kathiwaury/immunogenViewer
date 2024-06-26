% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotImmunogen.R
\name{plotImmunogen}
\alias{plotImmunogen}
\title{Plot protein features of one immunogen region}
\usage{
plotImmunogen(proteinDF, immunogen)
}
\arguments{
\item{proteinDF}{Protein DataFrame created by call to getProteinFeatures()}

\item{immunogen}{String, identifier name of immunogen}
}
\value{
A ggplot object
}
\description{
`plotImmunogen()` creates multiple ggplot objects within one figure. An Immunogen DataFrame is created by
filtering the Protein DataFrame for the relevant immunogen segment. A plot is created for each feature with annotations in
the Immunogen DataFrame. The amino acid sequence of the immunogen is shown on the x axis.
}
\examples{
proteinDF <- getProteinFeatures("P55087")
proteinDF <- addImmunogen(proteinDF, start=10, end=30, name="A12")
plotImmunogen(proteinDF, "A12")
}
