# clusterSystemsScript.R
#
# Purpose:
# Version: 1.0
# Date: April 1 2019
# Author: Rachel Silverstein
# License: MIT
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

if (FALSE) {
  systemDB_name <- "SysDB"
}


# ====  PACKAGES  ==============================================================
# Load all required packages.

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

# Define functions or source external files


# ====  PROCESS  ===============================================================

myDB <- fetchData(systemDB_name)
rootSysIDs <- SyDBgetRootSysIDs(myDB)

sys_names <- names(rootSysIDs)

# initialize a list of the genes in the different systems
systems <- vector(mode="list", length=length(sys_names))

names(systems) <- sys_names

# populate the list with vectors of genes
for (i in seq_along(sys_names)) {
  systems[[i]] <- SyDBgetSysSymbols(myDB, sys_names[i])
}


GTRD <- fetchData("GTRDgeneTFs")
GEO <- fetchData("GEOprofiles")
STRING <- fetchData("S")

result <- clusterSystems(systems,
                          distances = c("transcription_factor"),
                          customDistanceFn = NULL,
                          dataSources = NULL,
                          combineMatrices = TRUE,
                          printVennDiagrams = TRUE)


print(result)



# ====  TESTS  =================================================================
if (TRUE) {
  # Enter your function tests here...

}


# [END]
