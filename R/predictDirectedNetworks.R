# predictDirectedNetworks.R

# Purpose: To augment a system with genetic and physical interaction data in order to hypothesize regulatory relationshipd between components.
# Version: 1.0
# Date: 2019-03-21
# Author: Nada Elnour
# License: MIT

# Input: An excel spreadsheet of the system parsed using the BCB420-2019-resources scripts.
# Output: Graph hypotheses of possible regulatry networks in the system
# Dependencies: tibble 2.0.1; biomaRt 2.38.0; xlsx 0.6.1; readxl 1.3.1; dplyr 0.8.0.1; ggplot2 3.1.0; biogridr 0.0.0.9000; visNetwork 2.0.5.

# ====  GLOBAL PACKAGES ========================================================

source('./R/fetchData.R', echo=TRUE)
source('./R/SyDButils.R', echo=TRUE)
# ==============================================================================
#' Filter system by physically-interacting components and return their genetic interactions
#'
#' @param sysName A string specifying the full path to the excel sheet containing the system's components
#' @param criterion A string, either "stringent" or "relaxed", specifying whether only GGI between physical interactors should be selected.
#' @return A dataframe of system components and either their GGI between physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
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
