# STRINGExplore.R

#### getEdgeInteractions ####
# looks up the edge (from, to), made of 2 HGNC symbol strings, in the protActions dataset,
# returning an HTML string of the related protein actions in the STRING data system subset
# Input Example:
# getEdgeInteractions("VAMP3", "VAMP7", actionSet)
# Output Example : an HTML formatted string of the following form
# "nondirectional actions: binding <br /> actions that VAMP3 performs on VAMP7 : none found
# <br /> actions that VAMP7 performs on VAMP3 : catalysis; reaction"
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
#' @param genes (str) A 5 character string representing the name of
#' the desired curated gene system user would like to build a
#' high-confidence interaction network for.

#' @return A labeled list("nodes", "edges", "network") containing the visNetwork
#' object and the edge and node dataframes resulting from the setup of the network.
#' May be useful for downstream analysis by the end user.
#'
#'
#' @author \href{https://orcid.org/0000-0003-4762-8797}{Gabriela Morgenshtern} (aut)
#'
#' @import readr
#' @import BiocManager
#' @import biomaRt
#' @import visNetwork
#' @importFrom stats setNames
#' @importFrom magrittr %>%
#' @import BiocCheck
#' @import igraph
#'
#'
#' @examples
#' \dontrun{
#' STRINGExplore("PHALY")
#'}
#' @export

# it's better to ensure this doesn't rely on internet to run. make the input a gene list,
# and leverage the exported fn externally (3 inputs: geneSys, STRINGactions, STRINGedges)

STRINGExplore <- function(genes) {
  #### 1. Prepare / Load data:####
  # alternate / future route: add option to create edges from
  # high correlation scored expression data (Pearson rho >80%) and network those
  STRINGactions <- fetchData("STRINGactions")
  STRINGedges <- fetchData("STRINGedges0.8")
  geneSys <- list()

  dat <- fetchComponents(genes)
  # dat <- SyDBgetSysSymbols(fetchData("SysDB"), genes)
  geneSys <- unlist(dat)
  names(geneSys) <- c()

  actionSet <- STRINGactions[STRINGactions$protein1 %in% geneSys &
                               STRINGactions$protein2 %in% geneSys,]
  sel <- (STRINGedges$protein1 %in% geneSys) & (STRINGedges$protein2 %in% geneSys)
  geneSysEdges <- STRINGedges[sel, c("protein1", "protein2")]

  #### 2. igraph setup (visNetwork backbone)####
  sXG <-  igraph::graph.data.frame(geneSysEdges, directed = F)
  sXGgraph <-  igraph::simplify(sXG)

  #### 3. Calculate betweeness centrality of graph based on STRING scores ####
  bC <-  igraph::centr_betw(sXG)
  nodeBetw <- bC$res
  nodeBetw <- round(log(nodeBetw + 1)) + 1
  igraph::V(sXGgraph)$btwndegree <- nodeBetw

  #### 4. visNetwork setup from igraph object ####
  # 4.1 Nodes dataframe setup
  nodes <- igraph::get.data.frame(sXGgraph, what="vertices")
  nodes <- data.frame(id = nodes$name, title = nodes$name,
                      group = nodes$btwndegree, betweeness = nodes$btwndegree)
  setNames(nodes, c("id", "title", "betweeness", "betweeness centrality"))
  nodes <- nodes[order(nodes$id, decreasing = F),]

  # 4.2 Edges dataframe setup
  edges <- igraph::get.data.frame(sXGgraph, what="edges")
  edges$interactions <- ""
  for(i in 1:nrow(edges)){
    edges[i,]$interactions <- getEdgeInteractions(edges[i,]$from, edges[i,]$to, actionSet)
  }
  edges$title = paste("<p>",edges$from, "--",edges$to, ":", edges$interactions, "</p>")

  #### 5. Plot visNetwork ####
  network <- visNetwork::visNetwork(nodes, edges, height = "500px", width = "100%") %>%
    visNetwork::visOptions(selectedBy = "betweeness", highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visNetwork::visPhysics(stabilization = FALSE)

  #### 6. Return a labeled list for data transparency and future analysis ####
  networkData <- list("nodes" = nodes, "edges" = edges, "network" = network)
  return(networkData)
}

# [END]


