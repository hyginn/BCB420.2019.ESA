# ESA.R
#
# Purpose: Exploratory System Analysis tool project for BCB420 2019
# Version: 1.0.0
# Version history: 1.0.0
# Date: 25 March 2019
# Author: Deus Bajaj (deus.bajaj@mail.utoronto.ca)
# License: MIT
#
# Input: A biological system and a gene the user picks from the system
# Output: A histogram for the node degree distribution of the system gene network
# and a network plot for a gene from the biological system
# Dependencies: devtools, iGraph
#

# ====  FUNCTIONS  =============================================================

#'
#' \code{systemAnalysis()} Investigates genes from a biological system
#' and produces plots with gene network data
#'
#' @param gene A gene from a Biological system
#' @param sys A Biological system
#' @return a histogram for the node degree distribution of the system gene network
#' and a network plot for a gene from the biological system
#'
#' @author {Deus Bajaj} (aut)
#'
#' @examples
#' \dontrun{
#' # Retrieve a histogram for the node degrees of the PHALY gene network and
#' # get the network plot of thr AMBRA1 gene from the PHALY system
#' systemAnalysis("AMBRA1", "PHALY")
#' }
#'
#' @export

systemAnalysis <- function(gene, sys) {

  # Import STRINGedges from git repo at
  # https://github.com/deusbajaj/BCB420.2019.ESA
  # created from Dr. Steipe's BCB420.2019.ESA package
  # (forked directory, master repo can be accessed at
  # https://github.com/hyginn/BCB420.2019.ESA)

  # Load STRINGedges
  STRINGedges <- fetchData("STRINGedges0.9")

  myDB <- fetchData("SysDB")
  set <- SyDBgetSysSymbols(myDB, sys)

  # functionality adapted from Dr. Steipe's BCB420.2019.STRING package
  # (https://github.com/hyginn/BCB420.2019.STRING)
  # For this function we select edges for which both nodes are part of the set
  sel <- (STRINGedges$a %in% set) & (STRINGedges$b %in% set)
  xSetEdges <- STRINGedges[sel, c("a", "b")]

  sXGene <- igraph::graph_from_edgelist(matrix(c(xSetEdges$a,xSetEdges$b), ncol = 2))
  # degree distribution
  dg <- igraph::degree(sXGene)
  #Getting a histogram for the node degree distribution of the gene network
  x <- hist(dg, main = "Node degrees of the gene network",
       xlab = "Degree", ylab = "Counts")

  plot(x, col = "#A5CCF5", main = paste("Node degrees of the gene network"),
       sub = NULL, xlab = "Degree", ylab = "Counts")


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
