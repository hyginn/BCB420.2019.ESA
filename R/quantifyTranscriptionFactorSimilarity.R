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

# ==== FUNCTIONS ===============================================================

#' getSysGenes
#' Helper Function: Read system data from excel
#' @author \href{https://orcid.org/ORCID:0000-0002-1387-487X}{Rachel Woo} (aut)
#' @param filepath {string} to system excel sheet
#' Helper function to read in system data from excel
#' @return sysHGNC {list} A list of genes in the system
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

#' getSysTFdata
#' Helper Function: get TF data only for genes in system
#'
#' @param sysHGNC {list of string} list of genes from system
#' Accesses loaded TF data and filters for genes present in system
#' @return sysTF {list} A smaller TF geneList with only genes in system.
#' @author \href{https://orcid.org/ORCID:0000-0002-1387-487X}{Rachel Woo} (aut)
#' @examples
#' \donttest{
#' getSysTFdata(sysHGNC)
#' }
getSysTFdata <- function(sysHGNC) {
   # From Steipe, 2019 ESA readMe
  geneList <- fetchData("GTRDgeneTFs")
  # Only take elements from loaded data relevant to system
  sysTF <- geneList[sysHGNC]
  sysTF <- sysTF[lapply(sysTF,length)>0] #remove nulls

  return(sysTF)
}


#' makeDistanceMatrix
#' Helper function to make binary presence absence matrix
#'
#' @param sysTF {dataframe} output of getsysTFdata()
#' Makes a binary presence absence matrix which can be transformed into a
#' distance matrix.
#' @return matrix {dataframe} A matrix with rownames as TF and column names as system genes
#' In the matrix 1 means TF is associated with given gene. 0 means unassociated
#' @author \href{https://orcid.org/ORCID:0000-0002-1387-487X}{Rachel Woo} (aut)
#' @examples
#' \donttest{
#' makeDistanceMatrix(list)
#' }
makeDistanceMatrix <- function(sysTF) {
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

#' visualizeDendrogram
#' Helper Visualization Function
#'
#' @param matrix {dataframe} in the form of output of makeDistanceMatrix(). Binary matrix
#' Function which makes distanceMatrix, cluster and exports dendrogram
#' @return  {null} writes dendrogram in plots
#' @author \href{https://orcid.org/ORCID:0000-0002-1387-487X}{Rachel Woo} (aut)

#' @examples
#' \donttest{
#' visualizeDendrogram(matrix)
#' }
visualizeDendrogram <- function(matrix) {
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

#' quantifyTFSimilarity
#' Main controller function for this file. Effects all other functions.
#' @param filepath {string} to system excel file
#' Makes distance matrix, and cluster of shared TF presence for proteins in
#' the filepath denoted by system
#' @return {null} writes dendrogram of TF to plots
#' @export
#' @author \href{https://orcid.org/ORCID:0000-0002-1387-487X}{Rachel Woo} (aut)
#' @examples
#' \donttest{
#' quantifyTFSimilarity("BCB420-2019-System-PHALY-0.3.xlsx")
#' }
quantifyTFSimilarity <- function(filepath) {
  HGNC <- getSysGenes(filepath)
  sysTF <- getSysTFdata(HGNC)
  notInSysTF <- HGNC[!HGNC %in% names(sysTF)]
  matrix <- makeDistanceMatrix(sysTF)
  visualizeDendrogram(matrix)
}


