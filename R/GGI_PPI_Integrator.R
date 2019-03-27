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
#' @param filename A string specifying the full path to the excel sheet containing the system's components
#' @param mart A biomaRt mart to be queried for ENSEMBL ID conversion
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
#' @importFrom stats complete.cases
#' @export
getSysInteractions <-
  function(filename,
           mart = myMart,
           criterion = "stringent") {
    geneComp <-
      read_excel(filename, sheet = "component", skip = 1)
    # get my system's genes and their ENSEMBL IDs
    geneSym <- geneComp$sym[geneComp$molType != "concept"]

    att <- "ensembl_peptide_id"

    myGenes <- getBM(
      attributes = att,
      mart = myMart,
      filters = "hgnc_symbol",
      values = geneSym
    )
    myGenes <- myGenes[(myGenes$ensembl_peptide_id != ""), ]

    mySys <- convertToHGNC(myGenes, myMart)
    mySys <- getGeneticInteractome(mySys, criterion)

    return(mySys)
  }

#' Convert a system's components to HGNC-annotated physical interactors
#'
#' @param myGenes A list of system's ENSEMBL peptide IDs
#' @param mart A biomaRt mart to be queried for ENSEMBL ID conversion
#' @return The dataframe of physically interacting genes in \code{myGenes} according to \code{mart}
#' @import xlsx
#' @import readxl
#' @importFrom dplyr filter
#' @import biomaRt
#' @import ggplot2
#' @import biogridr
#' @import utils
#' @importFrom stats complete.cases
#' @export
convertToHGNC <- function(myGenes, mart) {

  interactions <-
    read.delim("../data/9606.protein.links.v11.0.txt", sep = " ")
  interactions$protein1 <-
    gsub("9606.", "", interactions$protein1)
  interactions$protein2 <-
    gsub("9606.", "", interactions$protein2)

  interactions <-
    interactions[(interactions$protein1 %in% myGenes &
                    interactions$protein2 %in% myGenes), ]
  interactions <- unique(interactions)

  interactions$protein1 <-
    recoverIDs(interactions$protein1, mart = mart)
  interactions$protein2 <-
    recoverIDs(interactions$protein2, mart = mart)

  ppi <-
    data.frame(interactions$protein1$sym, interactions$protein2$sym)
  ppi <- ppi[complete.cases(ppi), ]
  return(ppi)
}

#' Return BioGrid genetic interaction tags of physical interactors of the system
#'
#' @param mySys The dataframe output of \code{convertToHGNC}
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
getGeneticInteractome <- function(mySys, criterion) {
  myGenes <-
    toString(unique(c(
      as.character(mySys$interactions.protein1.sym),
      as.character(mySys$interactions.protein2.sym)
    )))

  myGenes <- gsub(", ", "|", myGenes)
  humInt <-
    bg("interactions") %>% bg_constrain(geneList = myGenes) %>% bg_get_results()
  humInt <-
    humInt[(
      humInt$experimental_system_type == "genetic" &
        humInt$organism_id_for_interactor_b == 9606
    ),]

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
             (
               match(gene1, mySys$interactions.protein1.sym) &
                 match(gene2, mySys$interactions.protein2.sym)
             ))
  } else if (criterion == "relaxed") {
    ggi <-
      filter(humInt,
             (
               match(gene1, mySys$interactions.protein1.sym) |
                 match(gene2, mySys$interactions.protein2.sym)
             ))
  }

  ggi <- ggi[complete.cases(ggi), ]
  return(ggi)
}

#' Generate the genetic interpretation map of BioGrid GGI tags of system's physical interactions.
#'
#' @return The dataframe mapping BioGrid GGI tag to its interpretation assuming that the system's components also interact physically.
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
#' @import visNetwork
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
#' myKey <- getKey("Tophie McGophie", "tmc@@hammertime.com", "SheTouchedThis")
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

# The following function is copied from Dr. Boris Steipe's STRING tool
# recoverIDs.R
#
# Purpose: To annnotate ENSEMBL protein IDs with gene names through BioMart
# Source: https://github.com/hyginn/BCB420.2019.STRING
# Version: 1.0
# Date: 2019-01-24
# Author: Boris Steipe
# License: MIT
# ==============================================================================
# NO SIDE EFFECTS:
# This script can be safely source()'d to define the functions it contains and
# install.packages()/run library() as required.
# All other code will not be executed unless this is done interactively.
# ==============================================================================

#' Try to recover IDs for ensp to sym mapping from biomart
#' @param  ensp (character) a vector of ensemble peptide IDs
#' @param mart (Mart)      an ensemble mart object
#' "HGNC" must exist in the global namespace
#' @return a dataframe with columns "ensp" containing the ensemble
#' peptide IDs of the input that could be mapped, and "sym",
#'  which contains the corresponding HGNC symbols, and rownames ensp.
#'
#' @export
recoverIDs <- function(ensp, mart) {
  # Purpose:
  #     Try to recover IDs for ensp to sym mapping from biomart
  # Parameters:
  #     ensp: (character) a vector of ensemble peptide IDs
  #     mart: (Mart)      an ensemble mart object
  #     "HGNC" must exist in the global namespace
  # Value:
  #     result: a dataframe with columns "ensp" containing the ensemble
  #             peptide IDs of the input that could be mapped, and "sym",
  #             which contains the corresponding HGNC symbols, and rownames
  #             ensp.

  # Note: to figure out the correct filters and attributes to use in
  #       querying a biomart, first fetch the filters and attributes with
  #       code like:
  #         filt <- biomaRt::listFilters(myMart)
  #         attr <- biomaRt::listAttributes(myMart)
  #       ... and then query:
  #         attr$name[grep("RefSeq", attr$description, ignore.case = TRUE)]

  # Define which attributes we want to fetch from biomart, and which columns
  # those match to in "HGNC":
  myAtt <- data.frame(biomart = c("uniprotswissprot",
                                  "refseq_mrna",
                                  "ucsc"),
                      HGNC =    c("UniProtID",
                                  "RefSeqID",
                                  "UCSCID"),
                      stringsAsFactors = FALSE)

  # Send off biomart query
  bm <- biomaRt::getBM(filters    =   "ensembl_peptide_id",
                       attributes = c("ensembl_peptide_id",
                                      myAtt$biomart),
                       values     = ensp,
                       mart       = mart)

  if (nrow(bm) > 0) {                   # at least one match was returned
    bm$sym <- rep(NA, nrow(bm))         # add a column to hold map results
    for (iCol in seq_len(ncol(bm))) {   # replace all "" with NA
      # Careful: combining logical vectors that can include NA is tricky.
      # Select elements that are neither already NA nor not-empty
      sel <- ( ! is.na(bm[ , iCol])) & (bm[ , iCol] == "")
      bm[sel, iCol] <- NA               # replace
    }
  }

  for (iAtt in seq_len(nrow(myAtt))) { # iterate over all requested attributes
    thisBmAtt <- myAtt$biomart[iAtt]
    thisHuAtt <- myAtt$HGNC[iAtt]
    if ( ! all(is.na(bm[ , thisBmAtt]))) {
      # some IDs were returned
      IDs <- bm[ , thisBmAtt]
      # get the symbol for a match, NA otherwise
      syms <- HGNC$sym[match(IDs, HGNC[ , thisHuAtt], incomparables = NA)]
      sel <- ( ! is.na(syms))
      bm$sym[sel] <- syms[sel] # Overwrite those that are not NA.
      # If there are multiple IDs returned for one row
      # in effect we return the last one that was
      # matched.
    }
  }
  # Post-process. Careful: ensemble_peptide_ids are not necessarily
  # unique in biomart output.
  bm <- bm[! is.na(bm$sym), c("ensembl_peptide_id", "sym")]
  bm <- bm[! duplicated(bm$ensembl_peptide_id), ]
  matchedIDs <- match(ensp, bm$ensembl_peptide_id)

  esMap <- data.frame(ensp = ensp,
                      sym = bm$sym[matchedIDs],
                      stringsAsFactors = FALSE)
  #rownames(esMap) <- esMap$ensp

  # drop NAs
  #esMap <- esMap[ ! is.na(esMap$sym), ]

  return(esMap)
}

# [END]

