# ESA.R
#
# Purpose: Exploratory System Analysis tool project for BCB420 2019
# Version: 1.0.0
# Version history: 1.0.0
# Date: 25 March 2019
# Author: Deus Bajaj (deus.bajaj@mail.utoronto.ca)
# License: MIT
#
# Input: The biological system PHALY and a gene the user picks from the system
# Output: A histogram for the node degree distribution of the PHALY gene network
# and a network plot for a gene from the PHALY system
# Dependencies: devtools, iGraph
#

# ====  FUNCTIONS  =============================================================

# function stub taken from Dr. Steipe's BCB420.2019.ESA package (https://github.com/hyginn/BCB420.2019.ESA)
fetchComponents <- function(sys) {
  # returns a fixed set of symbols.
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
  } else {
    s <- ""
  }
  return(s)
}

#'
#' \code{systemAnalysis()} Investigates genes from the system PHALY
#' and produces plots with gene network data
#'
#' @param gene A gene from the Biological system PHALY
#' @return a histogram for the node degree distribution of the PHALY gene network
#' and a network plot for a gene from the PHALY system
#'
#' @author {Deus Bajaj} (aut)
#'
#' @examples
#' systemAnalysis("AMBRA1")
#'
#' @export

systemAnalysis <- function(gene) {

  # Import STRINGedges from git repo at
  # https://github.com/deusbajaj/BCB420.2019.ESA
  # created from Dr. Steipe's BCB420.2019.ESA package
  # (forked directory, master repo can be accessed at
  # https://github.com/hyginn/BCB420.2019.ESA)

  # Load STRINGedges
  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "STRINGedges-2019-03-14.RData")
  load(url(myURL)) # loads STRING edges object

  set <- fetchComponents("PHALY") # Components for system PHALY

  # functionality adapted from Dr. Steipe's BCB420.2019.STRING package
  # (https://github.com/hyginn/BCB420.2019.STRING)
  # For this function we select edges for which both nodes are part of the set
  sel <- (STRINGedges$a %in% set) & (STRINGedges$b %in% set)
  xSetEdges <- STRINGedges[sel, c("a", "b")]

  sXGene <- igraph::graph_from_edgelist(matrix(c(xSetEdges$a,xSetEdges$b), ncol = 2))
  # degree distribution
  dg <- igraph::degree(sXGene)
  #Getting a histogram for the node degree distribution of the gene network
  hist(dg, col="#A5CCF5", main = "Node degrees of the gene network",
       xlab = "Degree", ylab = "Counts")


  gNet = STRINGedges[STRINGedges$a == gene|STRINGedges$b == gene ,]

  # functionality adapted from Dr. Steipe's BCB420.2019.STRING package
  # (https://github.com/hyginn/BCB420.2019.STRING)
  sXG <- igraph::graph_from_edgelist(matrix(c(gNet$a,gNet$b), ncol = 2))
  plot(sXG,
       layout = igraph::layout_with_fr(sXG),
       vertex.color=heat.colors(max(igraph::degree(sXG)+1))[igraph::degree(sXG)+1],
       vertex.size = 5.5 + (1.2 * igraph::degree(sXG)),
       vertex.label.cex = 0.5 + (0.025 * igraph::degree(sXG)),
       edge.width = 2,
       vertex.label = igraph::V(sXG)$name,
       vertex.label.family = "sans",
       vertex.label.cex = 0.9)
}

# [END]
