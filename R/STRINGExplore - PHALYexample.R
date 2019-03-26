#### Versions ####
# v 2.0 : functionalized example of tool output and usage
# v 1.0 : complete code for creating interactive viz, in script form
# v 0.1 : static coexpression igraph, using coexpression data;
#         separate list of protein actions related to edges
####
# Worked example of STRINGExplore with PHALY system

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

# Internal functions
source("./dev/STRINGExplore.R")
load("./data/STRINGedges.RData") #STRINGedges
load("./data/STRINGactions.RData") #STRINGactions

# fetchComponents is used in place of simply the system name to
# allow flexibility of input to the user-- for example, the possiblity
# of exploring a subsystem or simply a list of known HGNC symbols of interest,
# rather than one of the prepared curations
# PHALYnetwork
PHALYnetwork <- STRINGExplore(fetchComponents("PHALY"))
head(PHALYnetwork$nodes) # The dataframe of node information. Call to view the network nodes
head(PHALYnetwork$edges) # The dataframe of edge information. Call to view the network edges and protein interaction info
PHALYnetwork$network # The visNetwork object. Call to plot the network in your RStudio
