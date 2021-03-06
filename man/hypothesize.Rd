% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypothesize.R
\name{hypothesize}
\alias{hypothesize}
\title{hypothesize}
\usage{
hypothesize(network, ppi_ggi = NULL)
}
\arguments{
\item{network}{(dataframe) 3-column dataframe of physically-interacting
system components with their genetic interactions. The first two columns
denote the interacting pair; the third is the type of genetic interaction.}

\item{ppi_ggi}{(dataframe) An optional 3-column dataframe to specify subset
of \code{mySys} for which both PPI and GGI data is available. The first two
columns denote the interacting pair; the third is the type of genetic
interaction.}
}
\value{
(NULL) The function plots the graph of system components.
}
\description{
hypothesize
}
\examples{
\dontrun{
# Plot the graph of SLIGR components under stringent and relaxed conditions
mySys <- getSysInteractions("SLIGR", criterion = "stringent")
mySys2 <- getSysInteractions("SLIGR", criterion = "relaxed")
hypothesize(mySys) # draws hypothesis graph(s)
hypothesize(mySys2, mySys) # draws hypothesis graph(s)
}
}
\author{
Nada Elnour, \email{nada.elnour@mail.utoronto.ca}

\code{hypothesize} generates an annotated hypothesis graph whose nodes are
components of \code{mySys} and edges are GGI interpretations in either
\code{EMAP} or \code{GMAP}
}
