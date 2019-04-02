# quantifyTranscriptionFactorSimilarity.R
#
# Purpose: To quantify similarity between TF co-occurence and absence
# Version: 1.0
# Date: 2019-03-25
# Author: Rachel Woo
# License: MIT
#
# Input: An excel spreadsheet of the system parsed using the BCB420-2019-resources scripts.
# Output: Dengrogram of TF  co-occurance written to
#         ~/inst/img/TFOccuranceDendrogram.jpg
# Dependencies: xlsx, tidyverse, readxl
# ==============================================================================

# SIDE EFFECTS:
# Displays dendrogram

# # ====  PACKAGES  ==============================================================
# # Get packages
# if (requireNamespace("xlsx", quietly=TRUE)) {
#   install.packages("tidyverse", repos = "http://cran.us.r-project.org")
# }
# if (requireNamespace("devtools", quietly=TRUE)) {
#   install.packages("devtools", repos = "http://cran.us.r-project.org")
#   install.packages("Biostrings", repos = "http://cran.us.r-project.org")
#   install.packages("Rtools", repos = "http://cran.us.r-project.org")
# }
#
#
# # Load packages
# require(readxl, quietly = TRUE)


# ==== FUNCTIONS ===============================================================


#' Helper Function: Read system data from excel
#'
#' @param filepath to system excel sheet
#' Helper function to read in system data from excel
#' @return sysHGNC. A list of genes in the system
#' @examples
#' \donttest{
#' getSysGenes("BCB420-2019-System-PHALY-0.3.xlsx")
#' }
getSysGenes <- function(filepath) {
  # 1. Receive input data from excel spreadsheet and extract gene names with package tools
  sysComponent <- readxl::read_excel(filepath, sheet = "component", skip = 1)
  sysHGNC <- sysComponent$sym[sysComponent$molType == "protein"]
  #Get rid of NA
  sysHGNC <- sysHGNC[!is.na(sysHGNC)]
  return(sysHGNC)
}

#' Helper Function: get TF data only for genes in system
#'
#' @param sysHGNC list of genes from system
#' Accesses loaded TF data and filters for genes present in system
#' @return sysTF. A smaller TF geneList with only genes in system.
#' @examples
#' \donttest{
#' getSysTFdata(sysHGNC)
#' }
getSysTFdata <- function(sysHGNC) {
   # From Steipe, 2019 ESA readMe
   myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "geneList-2019-03-13.RData")
  load(url(myURL))  # loads GTRD geneList object
  # Only take elements from loaded data relevant to system
  sysTF <- geneList[sysHGNC]
  sysTF <- sysTF[lapply(sysTF,length)>0] #remove nulls

  return(sysTF)
}


#' Helper function to make binary presence absence matrix
#'
#' @param sysTF output of getsysTFdata()
#' Makes a binary presence absence matrix which can be transformed into a
#' distance matrix.
#' @return matrix. A matrix with rownames as TF and column names as system genes
#' In the matrix 1 means TF is associated with given gene. 0 means unassociated
#' @examples
#' \donttest{
#' makeMatrix(list)
#' }
makeMatrix <- function(sysTF) {
  #Make matrix with row and column names
  allTF <- unique(unlist(sysTF))
  matrix <- data.frame(matrix(0, nrow = length(names(sysTF)) , ncol = length(allTF)))
  rownames(matrix) <- names(sysTF)
  colnames(matrix) <- allTF

  #Add data to empty matrix
  for(gene in names(sysTF)) {
    geneIndex <- grep(gene, rownames(matrix))
    tfs <- sysTF[[gene]]
    for(i in 1:length(tfs)) {
      tf <- tfs[i]
      tfIndex <- grep(tf, colnames(matrix))
      matrix[geneIndex, tfIndex] <- 1
    }
  }
  return(matrix)
}

#' Helper Visualization Function
#'
#' @param matrix in the form of output of makeMatrix()
#' Function which makes distanceMatrix, cluster and exports dendrogram
#' @return  dendrogram in plots
#' @examples
#' \donttest{
#' visualize(matrix)
#' }
visualize <- function(matrix) {
  # Make binary distance matrix and cluster
  distanceMatrixBinary <- stats::dist(matrix, method = "binary")
  clusterBinary <- stats::hclust(distanceMatrixBinary)

  # Export Plots
  dendropath <- paste0(getwd(), "/inst/img/TFOccuranceDendrogram.jpg")
  heatpath <- paste0(getwd(), "/inst/img/TFOccuranceHeatmap.jpg")
  # Dendrogram
  graphics::plot(clusterBinary, main="Dendrogram of TF Co-occurence", xlab="Gene",
       ylab="Distance", sub="")
}

#' Main controller function for this file. Effects all other functions.
#' @param filepath to system excel file
#' Makes distance matrix, and cluster of shared TF presence for proteins in
#' the filepath denoted by system
#' @return dendrogram of TF
#' @export
#' @examples
#' \donttest{
#' quantifyTFSimilarity("BCB420-2019-System-PHALY-0.3.xlsx")
#' }
quantifyTFSimilarity <- function(filepath) {
  HGNC <- getSysGenes(filepath)
  sysTF <- getSysTFdata(HGNC)
  notInSysTF <- HGNC[!HGNC %in% names(sysTF)]
  matrix <- makeMatrix(sysTF)
  visualize(matrix)
}

#BCB420.2019.ESA::quantifyTFSimilarity("C:/Users/rwoo/Desktop/BCB420-2019-System-PHALY-0.3.xlsx")
#devtools::document("BCB420.2019.ESA")
#devtools::check("path/to/package/pkgname")
#devtools::build("path/to/package/pkgname")
#install.packages("path/to/package/packagename_version.tar.gz")
