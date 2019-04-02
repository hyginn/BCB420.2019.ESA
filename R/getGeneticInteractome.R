# getGeneticInteractome.R
#' @title getGeneticInteractome
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' \code{getGeneticInteractome} augments systems with BioGRID genetic
#' interaction annotations
#'
#' @param mySys The dataframe of PPI between system's components
#' @param criterion A string, either "stringent" or "relaxed", specifying
#' whether only GGI between physical interactors should be selected.
#'
#' @return The dataframe of system components and either their GGI between
#' physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
#'
#' @export
getGeneticInteractome <- function(mySys, criterion){
  # load BioGRID fetched GGI
  file <- system.file("extdata", "humInt.txt",
                      package = "BCB420.2019.ESA",
                      mustWork = TRUE)
  humInt <- read.delim(file, stringsAsFactors = FALSE)
  humInt <- humInt[complete.cases(humInt), ]

  if (criterion == "stringent") { # stringent => GGI between physical
                                  # interactors of system components

    interctorsInSystemID <- (humInt$A %in% mySys$a & humInt$B %in% mySys$b)

  } else if (criterion == "relaxed") { # stringent => all GGI in which system's
                                       # physical interactors are involved
    interctorsInSystemID <- (humInt$A %in% mySys$a | humInt$B %in% mySys$B)
  }

  ggi <- humInt[interctorsInSystemID,]
  ggi <- ggi[complete.cases(ggi),]

  return(ggi)
}

# [END]
