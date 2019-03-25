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

# load data for analysis
load("./data/STRINGedges.RData")
load("./data/STRINGactions.RData")



#### Analysis ####
# create adjacency matrix
xSet <- fetchComponents("PHALY")
sel <- (STRINGedges$protein1 %in% xSet) & (STRINGedges$protein2 %in% xSet)
xSetEdges <- STRINGedges[sel, c("protein1", "protein2")]

sXG <- igraph::graph_from_edgelist(matrix(c(xSetEdges$protein1,
                                            xSetEdges$protein2),
                                          ncol = 2,
                                          byrow = FALSE),
                                   directed = FALSE)

# calculate betweeness centrality based on coexpression scores

# Plot the network
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

# add Option for popup

# add protein action data to popup

