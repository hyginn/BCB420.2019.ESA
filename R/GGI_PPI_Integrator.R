# predictDirectedNetworks.R

# Purpose: To augment a system with genetic and physical interaction data in order to hypothesize regulatory relationshipd between components.
# Version: 1.0
# Date: 2019-03-21
# Author: Nada Elnour
# License: MIT

# Input: An excel spreadsheet of the system parsed using the BCB420-2019-resources scripts.
# Output: Graph hypotheses of possible regulatry networks in the system
# Dependencies: tibble 2.0.1; biomaRt 2.38.0; xlsx 0.6.1; readxl 1.3.1; dplyr 0.8.0.1; ggplot2 3.1.0; biogridr 0.0.0.9000; visNetwork 2.0.5
# ==============================================================================

# SIDE EFFECTS:
# This script imports biogridr which uses the deprecated function xml2::xml_find_one().

# ====  PACKAGES  ==============================================================
if (requireNamespace("biogridr", quietly = TRUE)) {
  library(biogridr)
} else {
  install.packages("devtools")
  devtools::install_github("npjc/biogridr")
  library(biogridr)
}

# Load all required packages.
require(xlsx, quietly = TRUE)
require(readxl, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(biomaRt, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(biogridr, quietly = TRUE)
require(visNetwork, quietly = TRUE)

# ==============================================================================
#' Filter system by physically-interacting components and return their genetic interactions
#'
#' @param sysName A string specifying the full path to the excel sheet containing the system's components
#' @param key a unique string key
#' @param criterion A string, either "stringent" or "relaxed", specifying whether only GGI between physical interactors should be selected.
#' @return A dataframe of system components and either their GGI between physical interactors should be selected (iff \code{criterion} == "stringent")
#' or all GGI of physical interactors (iff \code{relaxed} == "stringent")
#' @import xlsx
#' @import readxl
#' @importFrom dplyr filter
#' @import biomaRt
#' @import ggplot2
#' @import biogridr
#' @import utils
#' @include SyDButils.R fetchComponents.R fetchData.R
#' @importFrom stats complete.cases
#' @examples
#' \dontrun{
#' myKey <- getKey("Nada Elnour", "nada.elnour@@mail.utoronto.ca", "SLIGRESA")
#' mySys <- getSysInteractions("SLIGR", key = myKey, criterion = "stringent")
#' }
#'
#' @export
getSysInteractions <-
  function(sysName,
           key,
           criterion = "stringent") {
    if (is.null(sysName)) {
      stop("System not provided.\n")
    } else if (is.null(key)){
      stop("Provide valid key for biogrid ratification")
    }

    myURL <-
      paste0(
        "http://steipe.biochemistry.utoronto.ca/abc/assets/",
        "STRINGedges-2019-03-14.RData"
      )
    load(url(myURL))

    myDB <- fetchData("SysDB")

    geneComp <- SyDBgetSysSymbols(myDB, sysName)[[sysName]]

    interactions <- STRINGedges[(STRINGedges$a %in% geneComp &
                                   STRINGedges$b %in% geneComp),]
    interactions <-
      unique(interactions[complete.cases(interactions),])

    mySys <- getGeneticInteractome(mySys = interactions, criterion = criterion, key = key)

    return(mySys)
  }

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
      as.character(mySys$gene1),
      as.character(mySys$gene2)
    )))

  myGenes <- gsub(", ", "|", myGenes)
  humInt <-
    bg(access_point = "interactions", key = key) %>% bg_constrain(geneList = myGenes) %>% bg_get_results(.request = TRUE)
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
             (match(gene1, mySys$gene1) &
                match(gene2, mySys$gene2)))
  } else if (criterion == "relaxed") {
    ggi <-
      filter(humInt,
             (match(gene1, mySys$gene1) |
                match(gene2, mySys$gene2)))
  }

  ggi <- ggi[complete.cases(ggi),]
  return(ggi)
}

#' Generate the genetic interpretation map of BioGrid GGI tags of system's physical interactions.
#'
#' @return The dataframe mapping BioGrid GGI tag to its interpretation assuming that the system's components also interact physically.
#' @examples
#' EMAP <- makeEMAP()
#'
#' @export
makeEMAP <- function() {
  geneticInteractions <- c(
    "Dosage Growth Defect",
    "Dosage Lethality",
    "Dosage Rescue",
    "Negative Genetic",
    "Phenotypic Enhancement",
    "Phenotypic Suppression",
    "Positive Genetic",
    "Synthetic Growth Defect",
    "Synthetic Haploinsufficiency",
    "Synthetic Lethality",
    "Synthetic Rescue"
  )

  effects <- c(
    "inhibited by",
    "inhibited by",
    "equivalent/activates",
    "cis-regulates",
    #negative genetic
    "trans-regulates",
    "cis-regulates",
    "trans-regulates",
    #positive genetic
    "equivalent/activates",
    "equivalent/activates",
    "equivalent/activates",
    "inhibited by"
  )

  notes <- c(
    "",
    "",
    "in a redundant pathway",
    "if protein 2 is inhibitor, protein 1 inhibits protein 2; protein 2 is an activator, protein 1 activates protein 2; potentially antagonistic",
    #negative genetic
    "if protein 2 is inhibitor, protein 1 activates protein 2; protein 2 is an activator, protein 1 inhibits protein 2",
    "if protein 2 is inhibitor, protein 1 inhibits protein 2; protein 2 is an activator, protein 1 activates protein 2",
    "potentially synergistic",
    #positive genetic
    "in a redundant pathway",
    "in a redundant pathway",
    "in a redundant pathway",
    ""
  )

  EMAP <- data.frame(geneticInt = geneticInteractions,
                     effect = effects,
                     notes = notes)

  return(EMAP)

}

#' Generate the genetic interpretation map of BioGrid GGI tags.
#'
#' @return The dataframe mapping BioGrid GGI tag to its interpretation.
#' @examples
#' GMAP <- makeGMAP()
#'
#' @export
makeGMAP <- function() {
  geneticInteractions <- c(
    "Dosage Growth Defect",
    "Dosage Lethality",
    "Dosage Rescue",
    "Negative Genetic",
    "Phenotypic Enhancement",
    "Phenotypic Suppression",
    "Positive Genetic",
    "Synthetic Growth Defect",
    "Synthetic Haploinsufficiency",
    "Synthetic Lethality",
    "Synthetic Rescue"
  )

  effects <- c(
    "negative-parallels",
    "negative-parallels",
    "positive-parallels",
    "synergizes with",
    #negative genetic
    "synergizes with",
    "antagonizes",
    "antagonizes",
    #positive genetic
    "synergizes with",
    "synergizes with",
    "synergizes with",
    "antagonizes"
  )

  GMAP <- data.frame(geneticInt = geneticInteractions,
                     effect = effects)

  return(GMAP)

}

#' Generates an annotated hypothesis graph whose nodes are components of \code{mySys} and edges are GGI interpretations in either \code{EMAP} or \code{GMAP}
#'
#' @param network A dataframe of physically-interacting system components with their genetic interactions and
#' @param ppi_ggi An optional dataframe to specify subset of \code{mySys} for which both PPI and GGI data is available
#' @import visNetwork
#' @import xlsx
#' @import readxl
#' @importFrom dplyr filter
#' @import biomaRt
#' @import ggplot2
#' @import biogridr
#' @import utils
#' @include SyDButils.R fetchComponents.R fetchData.R
#' @importFrom stats complete.cases
#' @examples
#' \dontrun{
#' name <- readline("Name?")
#' email <- readline("email?")
#' project <- readline("project?")
#' myKey <- getKey(name, email, project)
#' mySys <- getSysInteractions("SLIGR", criterion = "stringent")
#' mySys2 <- getSysInteractions("SLIGR", criterion = "relaxed")
#' hypothesize(mySys)
#' hypothesize(mySys2, mySys)
#' }
#'
#' @export
hypothesize <-
  function(network,
           ppi_ggi = NULL) {
    EMAP <- makeEMAP()
    GMAP <- makeGMAP()
    visualizeInteractions(network, EMAP, ppi_ggi, GMAP)
  }

#' A helper function of \code{hypothesize} that plots the hypothesis graph
#'
#' @inheritParams hypothesize
#' @param gmap The dataframe obtained from \code{makeGMAP}.
#' @param emap The dataframe obtained from \code{makeEMAP}.
#'
#' @export
visualizeInteractions <- function(network, emap, ppi_ggi, gmap) {
  if (is.null(ppi_ggi))
  {
    network$gene1 <- as.factor(network$gene1)
    network$gene2 <- as.factor(network$gene2)
    allgenes <-
      as.factor(unique(c(
        as.character(network$gene1),
        as.character(network$gene2)
      )))
    nodes <-
      data.frame(id = allgenes,
                 allgenes,
                 label = unique(c(
                   as.character(network$gene1),
                   as.character(network$gene2)
                 )),
                 shape = 'circle')
    edges <-
      data.frame(
        from = network$gene1,
        to = network$gene2,
        label = emap$effect[match(network$interactionType, emap$geneticInt)],
        arrows = "to"
      )
  } else {
    network$gene1 <- as.factor(network$gene1)
    network$gene2 <- as.factor(network$gene2)

    ppi_ggi$gene1 <- as.factor(ppi_ggi$gene1)
    ppi_ggi$gene2 <- as.factor(ppi_ggi$gene2)

    allgenes <-
      as.factor(unique(c(
        as.character(network$gene1),
        as.character(network$gene2)
      )))

    allgenes2 <-
      as.factor(unique(c(
        as.character(ppi_ggi$gene1),
        as.character(ppi_ggi$gene2)
      )))

    nodes <-
      data.frame(
        id = allgenes,
        allgenes,
        label = unique(c(
          as.character(network$gene1),
          as.character(network$gene2)
        )),
        group = ifelse(allgenes %in% allgenes2 , "ppi-ggi", "ggi"),
        shape = 'circle'
      )

    edgeLabels <- c()

    for (i in 1:(length(network$interactionType))) {
      sel <-
        (
          as.character(network$gene1[i]) %in% as.character(ppi_ggi$gene1) &
            as.character(network$gene2[i]) %in% as.character(ppi_ggi$gene2)
        )
      if (sel) {
        idx <-
          which(
            as.character(ppi_ggi$gene1) == as.character(network$gene1[i]) &
              as.character(ppi_ggi$gene2) == as.character(network$gene2[i])
          )
        edgeLabels <-
          c(edgeLabels, as.character(emap$effect[match(ppi_ggi$interactionType[idx], emap$geneticInt)]))
      } else {
        edgeLabels <-
          c(edgeLabels, as.character(gmap$effect[match(network$interactionType[i], gmap$geneticInt)]))
      }
    }

    edges <-
      data.frame(
        from = network$gene1,
        to = network$gene2,
        label = edgeLabels,
        arrows = "to"
      )

  }

  visNetwork(nodes, edges) %>% visGroups(groupname = "ppi-ggi", color = "salmon") %>% visOptions(selectedBy = "group") %>% visEdges(color = "darkgray")
}

#' Generates a unique key for BioGrid API access
#'
#' @param my.name A string of the user's full name separated by space(s)
#' @param my.email A string matching a valid user email
#' @param my.project A string specifying a short project name (no spaces)
#' @return The unique key for user and project for data retrieval from BioGrid
#' @import biogridr
#' @examples
#' name <- readline("Name?")
#' email <- readline("email?")
#' project <- readline("project?")
#' myKey <- getKey(name, email, project)
#'
#' @export
getKey <- function(my.name, my.email, my.project) {
  my.name <- unlist(strsplit(my.name, " "))

  options(warn = -1)
  myKey <-
    bg_get_key(my.name[1], my.name[length(my.name)], my.email, my.project)
  options(warn = 0)
  return(myKey)
}

# [END]
