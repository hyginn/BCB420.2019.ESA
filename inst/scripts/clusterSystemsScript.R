# clusterSystemsScript.R
#
# Purpose: Script that generates the example results
#          used in my project page.
# Version: 1.0
# Date: April 1 2019
# Author: Rachel Silverstein
# License: MIT
#

# ==============================================================================

# NO SIDE EFFECTS:
# This script can be safely source()'d to define the functions it contains and
# install.packages()/run library() as required.
# All other code will not be executed unless this is done interactively.



# ====  PACKAGES  ==============================================================

if (require("cluster")) {
  install.packages(cluster)
}
library(cluster)

if (require("igraph")) {
  install.packages(igraph)
}
library(igraph)
if (require("VennDiagram")) {
  install.packages(VennDiagram)
}
library(VennDiagram)


# ====  FUNCTIONS  =============================================================

source("./R/clusterSystems.R")


# ====  PROCESS  ===============================================================

if (FALSE) {

  myDB <- fetchData("SysDB")
  rootSysIDs <- SyDBgetRootSysIDs(myDB)

  sys_names <- names(rootSysIDs)

  # initialize a list of the genes in the different systems
  systems <- vector(mode="list", length=length(sys_names))

  names(systems) <- sys_names

  # populate the list with vectors of genes
  for (i in seq_along(sys_names)) {
    systems[[i]] <- SyDBgetSysSymbols(myDB, sys_names[i])
  }

  # test clustering by variables individually
  clusterSystems(systems,
                 distances = c("transcription_factor"))

  clusterSystems(systems,
                 distances = c("expression_profile"))

  clusterSystems(systems,
                 distances = c("network_jaccard"))

  clusterSystems(systems,
                 distances = c("network_distance"))

  # test some combinations
  clusterSystems(systems,
                 distances = c("transcription_factor",
                               "network_jaccard",
                               "network_distance",
                               "expression_profile"),
                 combineMatrices = 'sum')

  clusterSystems(systems,
                 distances = c("transcription_factor",
                               "network_jaccard",
                               "network_distance",
                               "expression_profile"),
                 combineMatrices = 'product')

  clusterSystems(systems,
                 distances = c("transcription_factor",
                               "network_jaccard",
                               "network_distance",
                               "expression_profile"),
                 combineMatrices = 'minimum')

  clusterSystems(systems,
                 distances = c("transcription_factor",
                               "network_jaccard",
                               "network_distance",
                               "expression_profile"),
                 combineMatrices = 'maximum')

}





# [END]
