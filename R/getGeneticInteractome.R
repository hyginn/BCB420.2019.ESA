# getGeneticInteractome.R

#' Return BioGrid genetic interaction tags of physical interactors of the system
#'
#' @param mySys The dataframe of PPI between system's components
#' @param key  a unique string key
#' @param criterion A string, either "stringent" or "relaxed", specifying whether only GGI between physical interactors should be selected.
#' @return The dataframe of system components and either their GGI between physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
#'
#' @export
getGeneticInteractome <- function(mySys, criterion) {

  humInt <- fetchData("BioGRID")
  humInt <- humInt[humInt$type == "genetic", ]

  humInt$A <- toupper(humInt$A)
  humInt$B <- toupper(humInt$B)

  if (criterion == "stringent") {

    interctorsInSystemID <- (match(gene1, mySys$a) & match(gene2, mySys$b))
    ggi <- filter(humInt, interactorsInSystemID)

  } else if (criterion == "relaxed") {

    interctorsInSystemID <- (match(gene1, mySys$a) | match(gene2, mySys$b))
    ggi <- filter(humInt,interactorsInSystemID)

  }

  ggi <- ggi[complete.cases(ggi),]
  return(ggi)
}

# [END]
