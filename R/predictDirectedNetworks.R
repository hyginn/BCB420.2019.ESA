# predictDirectedNetworks.R
#
# Purpose: To augment a system with genetic and physical interaction data in order to hypothesize regulatory relationshipd between components.
# Version: 1.0
# Date: 2019-03-21
# Author: Nada Elnour
# License: MIT
#
# Input: An excel spreadsheet of the system parsed using the BCB420-2019-resources scripts.
# Output: Graph hypotheses of possible regulatry networks in the system
# Dependencies: tibble 2.0.1; biomaRt 2.38.0; xlsx 0.6.1; readxl 1.3.1; dplyr 0.8.0.1; ggplot2 3.1.0; biogridr 0.0.0.9000; visNetwork 2.0.5
# ==============================================================================

# SIDE EFFECTS:
# This script imports biogridr which uses the deprecated function xml2::xml_find_one().

# ====  PACKAGES  ==============================================================
# Load all required packages.
require(xlsx, quietly = TRUE)
require(readxl, quietly = TRUE)
require(dplyr, quietly = TRUE)
require(biomaRt, quietly = TRUE)
require(ggplot2, quietly = TRUE)
require(biogridr, quietly = TRUE)

my.name <- readline(prompt = "Enter full name separated by space(s): ")
my.name <- unlist(strsplit(my.name, " "))
my.email <- readline(prompt = "Enter email:")

my.project <-
  readline(prompt = "Enter project name without spaces:")

myKey <- bg_get_key(my.name[1], my.name[length(my.name)], my.email, my.project)

# ==============================================================================
getSysInteractions <-
  function(filename,
           intType = "genetic",
           mart = myMart) {
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

    # get interactome of interest

    mySys <- convertToHGNC(myGenes, myMart)
    mySys <- getGeneticInteractome(mySys)

    return(mySys)

  }

convertToHGNC <- function(myGenes, mart) {
  source("./R/recoverIDs.R")

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

  ##############
  interactions$protein1 <-
    recoverIDs(interactions$protein1, mart = mart)
  interactions$protein2 <-
    recoverIDs(interactions$protein2, mart = mart)

  #ggplot(interactions, aes(combined_score)) + geom_histogram() + theme_classic()

  ppi <-
    data.frame(interactions$protein1$sym, interactions$protein2$sym)
  ppi <- ppi[complete.cases(ppi), ]
  return(ppi)
}

getGeneticInteractome <- function(mySys) {
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

  ggi <-
    filter(humInt,
           (
             match(gene1, mySys$interactions.protein1.sym) &
               match(gene2, mySys$interactions.protein2.sym)
           ))

  ggi <- ggi[complete.cases(ggi), ]
  return(ggi)
}

makeEMAP <- function() {
  ##make epistasis map
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

  # In case no PPI data:
  notes <- c(
    "",
    "",
    "in a redundant pathway",
    "if protein 2 is inhibitor, protein 1 inhibits protein 2; protein 2 is an activator, protein 1 activates protein 2; potentially antagonistic",
    #negative genetic
    "if protein 2 is inhibitor, protein 1 activates protein 2; protein 2 is an activator, protein 1 inhibits protein 2",
    "if protein 2 is inhibitor, protein 1 inhibits protein 2; protein 2 is an activator, protein 1 activates protein 2",
    "potentially synergistic;",
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

##########################################
## Hypothesis Networks
hypothesize <- function(mySys, from = "physical") {
  EMAP <- makeEMAP()

  if (from == "physical") {
    visualizeInteractions(mySys, EMAP)
  }
}

##########################################
## Network Visualization
visualizeInteractions <- function(network, emap) {
  require(visNetwork, quietly = TRUE)
  network$gene1 <- as.factor(network$gene1)
  network$gene2 <- as.factor(network$gene2)
  allgenes <-
    as.factor(unique(c(
      as.character(network$gene1), as.character(network$gene2)
    )))
  nodes <-
    data.frame(id = allgenes,
               allgenes,
               label = unique(c(
                 as.character(network$gene1), as.character(network$gene2)
               )),
               shape = 'circle')
  edges <-
    data.frame(
      from = network$gene1,
      to = network$gene2,
      label = emap$effect[match(network$interactionType, emap$geneticInt)],
      arrows = "to"
    )
  visNetwork(nodes, edges, width = "100%")
}

###########################################
## Initialization

if (FALSE) {
  load("../BCB420-2019-resources/HGNC.RData")

  filename <- "./data/SLIGR.xlsx"
  ensembl <- useMart(biomart = "ensembl")
  human <- searchDatasets(mart = ensembl, pattern = "hsapiens")
  myMart <- useMart("ensembl", human$dataset)
  #######################################################################
  #######################################################################

  mySys <- getSysInteractions(filename)

  hypothesize(mySys, from = "physical")
}
#####################
#####################
# TODO:
# quantiative epistasis tools

# documentation
# extended hypotheses tool where even if genes are both physically and genetically interacting, they would be included in the hypothesis per genetic interaction results.
