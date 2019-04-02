# EDAtool.R
# Purpose: Exploratory Data Analysis tool project for BCB420 2019
# Version: 2.0.0
# Version history: 1.0.0
# Date: 02 April 2019
# Author: Yuhan Zhang(yuhan.zhang@mail.utoronto.ca)
# License: 
#MIT
#
# Input: set of genes of interest,could be fetched from a 5-letter system code like PHALY
# Output: heatmap of co-expression correlation, heatmap of GO semantic similarity of each pair of genes
# in the set. A scatterplot of correlation vs similarity and a fitted line within confidence level.
# 
# Dependencies: devtools,ggplot,ggcorrplot, GOSemSim
#

# ====  FUNCTIONS  =============================================================




# function stub taken from Dr. Steipe's BCB420.2019.ESA package (https://github.com/hyginn/BCB420.2019.ESA)
corGenes <- function(A, B, prf) {
  # Calculate pearson correlation between gene expression
  # profiles A and B in prf identified by the gene symbol.
  # A and B can be either gene symbol or index.
  
  r <- cor(prf[A, ], prf[B, ], use = "pairwise.complete.obs")
  return(r)
}

# used for sample codes and tests
# function stub taken from Dr. Steipe's BCB420.2019.ESA package (https://github.com/hyginn/BCB420.2019.ESA)
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
# ====  PACKAGES  ==============================================================
# Load all required packages.
#
# Use non-standard libraries with  package::function() idiom if possible.

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("GOSemSim", version = "3.8")
BiocManager::install("org.Hs.eg.db", version = "3.8")
install.packages("ggcorrplot")
# packages for calculating GO similarity
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="MF")


#'
#' \code{myPlot()} 
#' Investigates relationship of expression correlation and semantic similarity of input genes
#' and produces plots and summary of model info 
#'
#' @param geneSet A set of genes of interest to investigate
#' @return A heatmap of co-expression correlation of each pair of genes in the geneSet,
#' A heatmap of semantic similarity of each pair of genes
#' A scatterplot of correlation vs Go similarity 
#' summary of linear model built by correlation vs Go similarity
#'
#' @author {Yuhan Zhang} (aut)
#'
#' @examples
#' geneSet <- c("AMBRA1","ATG14","ATP2A1","ATP2A2","ATP2A3")
#' myPlot(geneSet)
#'
#' @export
  myPlot <- function(geneSet) {
    # Purpose:
    #     Collecting gene correlations and GO semantic similarity and build relationship model
    #     for each pair of genes to prove the hypothesis
    # Parameters:
    #     A : a set of genes of interets
    # Value:
    #     result: summary of the model withe scatterplot.
    
    # code ...
    # Loading necessary data
    myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                    "GEO-QN-profile-2019-03-24.rds")
    myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data
    myURL <- paste0("https://github.com/hyginn/",
                    "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
    load(url(myURL))  # loads HGNC data frame
    hsGO2 <- GOSemSim::godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)
    GeneA <- vector(mode = "character", length = 0)
    GeneB <- vector(mode = "character", length = 0)
    CorGene <- vector(mode = "numeric", length = 0)
    scoGO <- vector(mode = "numeric", length = 0)
    
    
    len <- length(geneSet)
    for (i in seq_along(geneSet)){
      for (j in i : len){
        a <- geneSet[i]
        b <- geneSet[j]
        if (i != j){ # different gene
         
         cor <- as.numeric(corGenes(a,b,myQNXP))
         
         scoreGO <- GOSemSim::geneSim(a, b, semData=hsGO2, measure="Wang", combine="BMA")$geneSim
        
        }else { # same gene
          cor <- 1
          scoreGO <- 1
          
        }
        GeneA <- c(GeneA, as.character(a))
        GeneB <- c(GeneB,as.character(b))
        CorGene <- c(CorGene, cor)
        scoGO <- c(scoGO,as.numeric(scoreGO))
      }
    }
    
    mydata <- data.frame(GeneA, GeneB, CorGene,scoGO)
    corM <- matrix(nrow = len, ncol = len)
    goM <- matrix(nrow = len,ncol = len)
    
    colnames(corM) <- unique(mydata$GeneB)
    rownames(corM) <- unique(mydata$GeneA)
    colnames(goM) <- unique(mydata$GeneB)
    rownames(goM) <- unique(mydata$GeneA)

    for (i in 1:nrow(mydata)) {
      corM[
        mydata$GeneA[i],
        mydata$GeneB[i]
        ] <- CorGene[i]
      corM[
        mydata$GeneB[i],
        mydata$GeneA[i]
        ] <- CorGene[i]
      goM[
        mydata$GeneA[i],
        mydata$GeneB[i]
        ] <- scoGO[i]
      goM[
        mydata$GeneB[i],
        mydata$GeneA[i]
        ] <- scoGO[i]
    }
    corM[is.na(corM)] <- 0
    heatmap1 <- ggcorrplot::ggcorrplot(corM, hc.order = TRUE, type = "lower",
                           outline.col = "white", lab = TRUE)
    heatmap2 <- ggcorrplot::ggcorrplot(goM, hc.order = TRUE, type = "lower",
                           outline.col = "white", lab = TRUE)
    print(heatmap1)
    print(heatmap2)
      pairFrame <- mydata[which(mydata$GeneB != mydata$GeneA),]
      Correlation <- pairFrame$CorGene
      goSimilarity <- pairFrame$scoGO
      model <- lm(goSimilarity~Correlation)
      lmplot <- ggplot2::ggplot(pairFrame, ggplot2::aes(x=Correlation, y=goSimilarity)) + 
        ggplot2::geom_point()+
        ggplot2::geom_smooth(method=lm)+
        ggplot2::ggtitle("Co-expression Correlation VS GO Semantic Similarity") +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(name = "Semantic Similarity") +
        ggplot2::scale_y_continuous(name = "Co-expression Correlation")
     print(lmplot)
    print(summary(model))
      
    
    return(mydata)
  }
  
geneSet <- c("AMBRA1","ATG14","ATP2A1","ATP2A2","ATP2A3")
myPlot(geneSet)
