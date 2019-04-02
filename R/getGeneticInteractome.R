# getGeneticInteractome.R

#' Return BioGrid genetic interaction tags of physical interactors of the system
#'
#' @param mySys The dataframe of PPI between system's components
#' @param criterion A string, either "stringent" or "relaxed", specifying whether only GGI between physical interactors should be selected.
#' @return The dataframe of system components and either their GGI between physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
#'
#' @export
getGeneticInteractome <- function(mySys, criterion){
  # load BioGRID fetched GGI
  file <- system.file("extdata", "humInt.txt", package = "BCB420.2019.ESA", mustWork = TRUE)
  humInt <- read.delim(file)
  humInt <- humInt[complete.cases(humInt), ]

  if (criterion == "stringent") { # stringent => GGI between physical interactors of system components

    interctorsInSystemID <- (humInt$A %in% mySys$a & humInt$B %in% mySys$b)
    ggi <- filter(humInt, interctorsInSystemID)

  } else if (criterion == "relaxed") { # stringent => all GGI in which a system's physical interactors are involved

    interctorsInSystemID <- (humInt$A %in% mySys$a | humInt$B %in% mySys$B)
    ggi <- filter(humInt, interctorsInSystemID)

  }

  ggi <- ggi[complete.cases(ggi),]
  return(ggi)
}

# [END]
