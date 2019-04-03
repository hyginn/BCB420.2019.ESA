# predictDirectedNetworks
#
#' @title getSysInteractions
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' \code{getSysInteractions} filter system by physically-interacting components
#' and return their genetic interactions
#'
#' @param sysName (string|vector of strings) Specifies the systems of interest
#' @param criterion (string) Either "stringent" or "relaxed", specifying if
#' only GGI between physical interactors should be selected.
#'
#' @return (dataframe) A 3-column dataframe of system components and either
#' their GGI between physical interactors should be selected
#' (iff \code{criterion} == "stringent") or all GGI of physical interactors
#' (iff \code{relaxed} == "stringent"). The first two columns denote the
#' interacting pair; the third is the type of genetic interaction.
#'
#' @importFrom dplyr filter
#' @import utils
#' @include SyDButils.R fetchData.R
#' @importFrom stats complete.cases
#' @examples
#' # Get SLIGR physical interactors and return the genetic interactions between
#' # them
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
