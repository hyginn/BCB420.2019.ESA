# getGeneticInteractome.R

#' Return BioGrid genetic interaction tags of physical interactors of the system
#'
#' @param mySys The dataframe of PPI between system's components
#' @param key  a unique string key
#' @param criterion A string, either "stringent" or "relaxed", specifying whether only GGI between physical interactors should be selected.
#' @return The dataframe of system components and either their GGI between physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
#' @import xlsx
#' @import readxl
#' @importFrom dplyr filter
#' @import biomaRt
#' @import ggplot2
#' @import biogridr
#' @import utils
#' @importFrom stats complete.cases
#' @export
getGeneticInteractome <- function(mySys, criterion, key) {
  myGenes <-
    toString(unique(c(
      as.character(mySys$a),
      as.character(mySys$b)
    )))

  myGenes <- gsub(", ", "|", myGenes)
  humInt <-
    bg(access_point = "interactions", key = key) %>% bg_constrain(geneList = myGenes) %>% bg_get_results()
  humInt <-
    humInt[(
      humInt$experimental_system_type == "genetic" &
        humInt$organism_id_for_interactor_b == 9606
    ), ]

  humInt <-
    data.frame(
      humInt$official_symbol_for_interactor_a,
      humInt$official_symbol_for_interactor_b,
      humInt$experimental_system_name
    )

  colnames(humInt) <- c("gene1", "gene2", "interactionType")
  humInt$gene1 <- toupper(humInt$gene1)
  humInt$gene2 <- toupper(humInt$gene2)

  if (criterion == "stringent") {
    ggi <-
      filter(humInt,
             (match(gene1, mySys$a) &
                match(gene2, mySys$b)))
  } else if (criterion == "relaxed") {
    ggi <-
      filter(humInt,
             (match(gene1, mySys$a) |
                match(gene2, mySys$b)))
  }

  ggi <- ggi[complete.cases(ggi),]
  return(ggi)
}

# [END]
