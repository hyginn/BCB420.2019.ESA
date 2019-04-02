# predictDirectedNetworks.R

source('./R/fetchData.R', echo=TRUE)
source('./R/SyDButils.R', echo=TRUE)
# ==============================================================================
#' @title getSysInteractions
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' \code{getSysInteractions} filter system by physically-interacting components
#' and return their genetic interactions
#'
#' @param sysName A string specifying the full path to the excel sheet
#' containing the system's components
#' @param criterion A string, either "stringent" or "relaxed", specifying if
#' only GGI between physical interactors should be selected.
#'
#' @return A dataframe of system components and either their GGI between
#' physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
#'
#' @importFrom dplyr filter
#' @import utils
#' @include SyDButils.R fetchData.R
#' @importFrom stats complete.cases
#' @examples
#' mySys <- getSysInteractions("SLIGR", criterion = "stringent")
#'
#' @export
getSysInteractions <-function(sysName, criterion = "stringent") {
    # check if a system name is given
    if (is.null(sysName)) {
      stop("System not provided.\n")
    }

    # fetch the system's components
    STRINGedges <- as.data.frame(fetchData("STRINGedges0.9"))

    myDB <- fetchData("SysDB")

    systems <- SyDBgetSysSymbols(myDB, sysName)
    geneComp <- character()

    if (length(sysName) == 1) { # if only one system is to be fetched
      geneComp <- SyDBgetSysSymbols(myDB, sysName)[[sysName]]
    } else {
      for (system in sysName) {
        geneComp <- c(geneComp, systems[[system]])
      }
    }

    # map to STRING dataset to get PPI between systems' components
    interactions <- as.data.frame(STRINGedges[(STRINGedges$a %in% geneComp &
                                   STRINGedges$b %in% geneComp),])
    interactions <-
      unique(interactions[complete.cases(interactions),])

    # augment PPI with GGI data
    mySys <- getGeneticInteractome(mySys = interactions, criterion = criterion)

    return(mySys)
  }

# [END]
