# STRINGExplore.R
#### Helper functions ####
if (! requireNamespace("readr")) {
  install.packages("readr")
}
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

if (! requireNamespace("devtools")) {
  install.packages("devtools")
}

# Question for B Steipe: why didn't this download the mapping tool script into /inst?

if (! requireNamespace("devtools")) {
  install.packages("devtools")
  devtools::install_github("hyginn/BCB420.2019.STRING")
}

if (! requireNamespace("BiocCheck")) {
  BiocManager::install("BiocCheck")
}
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
  } else {
    s <- ""
  }
  return(s)
}

#### getEdgeInteractions
# looks up the edge (from, to), made of 2 HGNC symbol strings, in the protActions dataset,
# returning an HTML string of the related protein actions in the STRING data system subset
# Input Example:
# getEdgeInteractions("VAMP3", "VAMP7", actionSet)
# Output Example : an HTML formatted string of the following form
# "nondirectional actions: binding <br /> actions that VAMP3 performs on VAMP7 : none found
# <br /> actions that VAMP7 performs on VAMP3 : catalysis; reaction"
#####
getEdgeInteractions <- function(from, to, sysActions){

  interactionSet <- sysActions[((sysActions$protein1==from & sysActions$protein2==to) |
                                  (sysActions$protein1==to & sysActions$protein2==from)),]
  # # Create HTML
  ndActions <- paste(interactionSet[interactionSet$is_directional == FALSE,]$mode, collapse = "; ")
  if(ndActions == "") {
    ndActions <- "none found"
  }
  nondirectional <- paste("nondirectional actions:", paste(ndActions, collapse = " "))


  # #from-->to directional actions
  ftActions <- paste(interactionSet[interactionSet$is_directional == TRUE &
                                      interactionSet$a_is_acting == TRUE,]$mode, collapse = "; ")
  if(ftActions == "") {
    ftActions <- "none found"
  }
  fromTo <- paste("actions that", from, "performs on", to, ":", paste(ftActions, collapse = " "))

  #to-->from directional actions

  tfActions <- paste(interactionSet[interactionSet$is_directional == TRUE &
                                      interactionSet$a_is_acting == FALSE,]$mode, collapse = "; ")
  if(tfActions == "") {
    tfActions <- "none found"
  }
  toFrom <- paste("actions that", to, "performs on", from, ":", paste(tfActions, collapse = " "))

  interactions <- paste(nondirectional, "<br />" , fromTo, "<br />", toFrom)
  return(interactions)
}

#' STRINGExplore
#'
#' \code{STRINGExplore} Utilizes the combined score and protein action information
#'  available on STRING to create an interactive network of high confidence (>0.8) interaction
#'
#' @section Details: The STRING edges on their don't provide a good starting
#'  point for hypothesizing about the reasons for the existence of high-
#'  confidence interactions. By combining the easy-to-navigate, interactive
#'  network from the visNetwork package with tooltips of information from
#'  high-confidence (>0.8) data from STRING protein_actions dataset, STRINGExplore
#'  aims to ease the process of deciding on future research avenues. The graph layout takes
#'  into consideration the betweeness centrality of nodes, and so when inputting a gene list
#'  from a curated system it's highly likely that highly-centrality nodes will impart some
#'  important functionality to the system.
#'
#' @param geneList (char set) HGNC symbols representing the desired list of genes
#'  user would like to build a high-confidence interaction network for.
#' @param globalEdges (dataframe). The set of high-confidence edges available for
#' querying for geneList members.Each column is an HGNC symbol.
#' The expected characteristics can be found in R/mapSTRING.R for STRINGedges.

#' @param globalActions (dataframe) <description>. The set of high-confidence protein
#' actions available for querying for geneList members.Each column is an HGNC symbol.
#' The expected characteristics can be found in R/mapSTRING.R for STRINGactions.

#' @return A labeled list("nodes", "edges", "network") containing the visNetwork
#' object and the edge and node dataframes resulting from the setup of the network.
#' May be useful for downstream analysis by the end user.
#'
#'
#' @author \href{https://orcid.org/0000-0003-4762-8797}{Gabriela Morgenshtern} (aut)
#'
#'
#' @examples
#' # A user inputs a list of HGNC gene symbols and STRINGExplore will look them
#' up in the prepared high-confidence STRING datasets for edges and protein actions,
#' formatting an interactive network complete with protein interaction data
#' STRINGExplore(fetchComponents("PHALY"))
#'
#' @export

STRINGExplore <- function(geneList, globalEdges, globalActions) {
  #### 1. Prepare / Load data:####
  # alternate / future route: add option to create edges from
  # high correlation scored expression data (Pearson rho >80%) and network thos
  geneSys <- fetchComponents("PHALY")
  actionSet <- globalActions[globalActions$protein1 %in% geneSys &
                               globalActions$protein2 %in% geneSys,]
  sel <- (globalEdges$protein1 %in% geneSys) & (globalEdges$protein2 %in% geneSys)
  geneSysEdges <- globalEdges[sel, c("protein1", "protein2")]

  #### 2. Calculate betweeness centrality of graph based on STRING scores ####
  bC <- centr_betw(sXG)
  nodeBetw <- bC$res
  nodeBetw <- round(log(nodeBetw + 1)) + 1

  #### 3. igraph setup (visNetwork backbone)####
  sXG <- igraph::graph.data.frame(geneSysEdges, directed = F)
  sXGgraph <- igraph::simplify(sXG)
  V(sXGgraph)$btwndegree <- nodeBetw

  #### 4. visNetwork setup from igraph object ####
  # 4.1 Nodes dataframe setup
  nodes <- get.data.frame(sXGgraph, what="vertices")
  nodes <- data.frame(id = nodes$name, title = nodes$name,
                      group = nodes$btwndegree, betweeness = nodes$btwndegree)
  setNames(nodes, c("id", "title", "betweeness", "betweeness centrality"))
  nodes <- nodes[order(nodes$id, decreasing = F),]

  # 4.2 Edges dataframe setup
  edges <- get.data.frame(graph, what="edges")
  edges$interactions <- ""
  for(i in 1:nrow(edges)){
    edges[i,]$interactions <- getEdgeInteractions(edges[i,]$from, edges[i,]$to, actionSet)
  }
  edges$title = paste("<p>",edges$from, "--",edges$to, ":", edges$interactions, "</p>")

  #### 5. Plot visNetwork ####
  network <- visNetwork(nodes, edges, height = "500px", width = "100%") %>%
    visOptions(selectedBy = "betweeness", highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visPhysics(stabilization = FALSE)

  #### 6. Return a labeled list for data transparency and future analysis ####
  networkData <- list("nodes" = nodes, "edges" = edges, "network" = network)
  return(networkData)
}

# [END]
