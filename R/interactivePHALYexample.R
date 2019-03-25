#### Versions ####
# v 0.1 : static coexpression igraph, using coexpression data;
#         separate list of protein actions related to edges
####
# Worked example ofSTRINGCoExplore with PHALY system

#### Setup and Dependencies ####
# Package requirements
if (! require(igraph, quietly=TRUE)) {
  install.packages("igraph")
  library(igraph)
}

if (!require(RColorBrewer, quietly=TRUE)) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

if (!require(magrittr, quietly=TRUE)) {
  install.packages("magrittr")
  library(magrittr)
}

if (!require(visNetwork, quietly=TRUE)) {
  install.packages("visNetwork")
  library(visNetwork)
}

# Helper functions
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

# looks up the edge (from, to) in the protActions dataset, returning
# an HTML string of the related protein actions in the STRING data system subset
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
getEdgeInteractions("VAMP3", "VAMP7", actionSet)

# load data for analysis
# alternate: expression data, find the edges w coexpression >60%, network those
load("./data/STRINGedges.RData")

load("./data/STRINGactions.RData")

#### Load the expression profiles:####
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "GEO-QN-profile-2019-03-24.rds")
myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data

corGenes <- function(A, B, prf) {
  # Calculate pearson correlation between gene expression
  # profiles A and B in prf identified by the gene symbol.
  # A and B can be either gene symbol or index.

  r <- cor(prf[A, ], prf[B, ], use = "pairwise.complete.obs")
  return(r)
}

# Returns a matrix of coexpression scores for the given HGNC symbol list
# and associated quantile-normalized expression dataset
coexpressionMatrix <- function(geneList, prf) {
  # TO DO: colnames, also, data frame is probably easiest here
  coexMatrix <- matrix(0, nrow = length(testList),
                       ncol = length(testList), dimnames = list(testList))
  for(gene in geneList){
  # TO DO: we want to add the results from this function to the coexpression matrix
    lapply(geneList[geneList != gene], corGenes, A = gene, prf = prf)
  }
}

#### Analysis ####
# create adjacency matrix
xSet <- fetchComponents("PHALY")
actionSet <- STRINGactions[STRINGactions$protein1 %in% xSet &
                             STRINGactions$protein2 %in% xSet,]
sel <- (STRINGedges$protein1 %in% xSet) & (STRINGedges$protein2 %in% xSet)
xSetEdges <- STRINGedges[sel, c("protein1", "protein2")]

# sXG <- igraph::graph_from_edgelist(matrix(c(xSetEdges$protein1,
#                                             xSetEdges$protein2),
#                                           ncol = 2,
#                                           byrow = FALSE),
#                                    directed = FALSE)

# calculate betweeness centrality based on STRING/coexpression scores
bC <- centr_betw(sXG)
nodeBetw <- bC$res
nodeBetw <- round(log(nodeBetw + 1)) + 1

sXG <- igraph::graph.data.frame(xSetEdges, directed = F)
sXGgraph <- igraph::simplify(sXG)
V(sXGgraph)$btwndegree <- nodeBetw


#### Plot the network (igraph - from BCH441)####
oPar <- par(mar= rep(0,4)) # Turn margins off
set.seed(112358)
plot(sXG,
     layout = igraph::layout_with_fr(sXG),
     vertex.color=heat.colors(max(igraph::degree(sXG)+1))[igraph::degree(sXG)+1],
     vertex.size = 1.5 + (1.2 * igraph::degree(sXG)),
     vertex.label.cex = 0.2 + (0.025 * igraph::degree(sXG)),
     edge.width = 2,
     vertex.label = igraph::V(sXG)$name,
     vertex.label.family = "sans",
     vertex.label.cex = 0.9)
set.seed(NULL)
par(oPar)


#### visNetwork plot (from visNetwork code examples) ####
# Nodes dataframe setup
nodes <- get.data.frame(sXGgraph, what="vertices")
nodes <- data.frame(id = nodes$name, title = nodes$name,
                    group = nodes$btwndegree, betweeness = nodes$btwndegree)
setNames(nodes, c("id", "title", "betweeness", "betweeness centrality"))
nodes <- nodes[order(nodes$id, decreasing = F),]

#### Edges dataframe setup ####
edges <- get.data.frame(graph, what="edges")
edges$interactions <- ""
for(i in 1:nrow(edges)){
  edges[i,]$interactions <- getEdgeInteractions(edges[i,]$from, edges[i,]$to, actionSet)
}
edges$title = paste("<p>",edges$from, "--",edges$to, ":", edges$interactions, "</p>")

#### Plot visNetwork ####
visNetwork(nodes, edges, height = "500px", width = "100%") %>%
  visOptions(selectedBy = "betweeness", highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visPhysics(stabilization = FALSE)



