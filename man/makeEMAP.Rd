% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeMaps.R
\name{makeEMAP}
\alias{makeEMAP}
\title{Generate the genetic interpretation map of BioGrid GGI tags of system's
 physical interactions.}
\usage{
makeEMAP()
}
\value{
(dataframe) An 11-by-3 dataframe mapping BioGrid GGI tag to its
 interpretation assuming that the system's components also interact
 physically. The first column contains the official BioGRID GGI tags; the
 second contains interpreted relationships under the assumption; the third
  contains notes on interpretation.
}
\description{
Generate the genetic interpretation map of BioGrid GGI tags of system's
 physical interactions.
}
\examples{
\dontrun{
EMAP <- makeEMAP() # generates and loads the EMAP
}
}
