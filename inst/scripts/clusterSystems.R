# scriptTemplate.R
#
# Purpose:
# Version:
# Version history:
# Date:
# Author:
# License:
#
# Input:
# Output:
# Dependencies:
#
# ToDo:
# Notes:
#
# ==============================================================================

# NO SIDE EFFECTS:
# This script can be safely source()'d to define the functions it contains and
# install.packages()/run library() as required.
# All other code will not be executed unless this is done interactively.


# ====  PARAMETERS  ============================================================
#
#      This script uses guard-blocks that prevent execution of
#      code that should not be executed when the entire script
#      is sourced. Thus it can be source()'d to load its functions,
#      or executed interactively.
#
if (TRUE) {

  # list the names of the systems to be used in analysis
  sys_names <- c("PHALY", "ROUSI")

  # load data frame of HGNC symbols, "HGNC"
  myURL <- paste0("https://github.com/hyginn/",
                  "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))  # loads HGNC data frame

}


# ====  PACKAGES  ==============================================================
# Load all required packages.
#
# Use non-standard libraries with  package::function() idiom if possible.

if (!require(cluster, quietly=TRUE)) {
  install.packages("cluster")
}
library(cluster)

if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
}

if (!require(devtools)) {
  install.packages("devtools")
}
library(devtools)

if (!require(igraph)) {
  install.packages("igraph")
}
library(igraph)


# ====  FUNCTIONS  =============================================================

# Define functions or source external files
if (TRUE) {

  # load pairwise similarity functions
  source("./R/expr_sim.R")
  source("./R/tf_sim.R")
  source("./R/net_sim.R")
  # load function to make similarity matrices
  source("./R/make_matrix.R")
}

# remove before pull request
fetchComponents <- function(sys) {
  # returns a fixed set of symbols.
  # Function stub for development purposes only.
  if (sys == "PHALY") {
    s <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
           "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
           "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
           "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
           "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
           "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
           "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
           "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
           "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
           "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
           "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
           "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
           "VPS41", "VTI1B", "YKT6")
  } else if (sys == "ROUSI") {
    s <- c("ROBO3", "SLIT2", "ROBO2", "ROBO1", "SLIT3", "MYO9B", "RHOA")
  } else {
    s <- ""
  }
  return(s)
}


# ====  PROCESS  ===============================================================

if (TRUE) {

  # === Fetch the System Components ======

  # initialize a list of the different systems
  systems <- vector(mode="list", length=length(sys_names))
  names(systems) <- sys_names

  # populate the list with vectors of genes
  for (i in seq_along(sys_names)) {
    systems[[i]] <- fetchComponents(sys_names[i])
  }

  # make a single vector containing all of the genes from all of the systems together
  all_genes <- unlist(systems)


  # ===== Load in the data by which we are going to cluster the genes =======

  # 1. EXPRESSION PROFILES
  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "GEO-QN-profile-2019-03-24.rds")
  myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data


  # 2. TRANSCRIPTION FACTORS
  # list of TFs that have a binding site upstream of gene
  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "geneList-2019-03-13.RData")
  load(url(myURL))  # loads GTRD geneList object

  # 3. STRING EDGES
  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "STRINGedges-2019-03-14.RData")
  load(url(myURL))  # loads STRING edges object
  # convert the STRINGedges object into an igraph object
  STRINGgraph <- graph_from_edgelist(as.matrix(STRINGedges[,1:2]))
  # The string graph does not contain any unconnected vertices, only interactions,
  # so any genes that have no annotated interactions will be missing.
  # Add these missing genes as unconnected nodes in the graph.
  string_genes <- unique(c(STRINGedges$a, STRINGedges$b))
  missing_genes <- setdiff(all_genes, string_genes)
  STRINGgraph <- add.vertices(STRINGgraph, nv = length(missing_genes), name = missing_genes)


  # ======= Create Similarity Matrices =======

  # pass the corresponding similarity function for each matrix
  tf_sim_matrix <- make_matrix(all_genes, tf_sim) # transcription factor matrix
  expr_sim_matrix <- make_matrix(all_genes, expr_sim) # expression matrix

  # convert the similarity matrices to distance (dissimilarity) matrices
  # or just directly make distance matrix for the string data
  tf_dist_matrix <- as.dist(sqrt(1 - tf_sim_matrix))
  expr_dist_matrix <- as.dist(sqrt(1 - expr_sim_matrix))


  string_matrix <- igraph::distances(STRINGgraph, v = all_genes, to = all_genes)
  # remove all the infinite values generated by unconnected nodes
  # get a number one larger than the largest finite distance in the matrix
  large_num <- max(string_matrix[is.finite(string_matrix)]) + 1
  # replace any infinite distances with this distance
  sel <- is.infinite(string_matrix)
  string_matrix[sel] <- large_num

  string_dist_matrix <- as.dist(string_matrix)



  # ====== Cluster ========














}

# ====  TESTS  =================================================================
if (TRUE) {
  # Enter your function tests here...

}


# [END]
