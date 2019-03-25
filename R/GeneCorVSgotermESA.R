# ESATool.R
#
# Purpose: To show whether gene expression similarity contain linear relationship with related GO terms, 
# Version:1.0
# Version history:
# Date:
# Author:Yuhan Zhang
# License:
#
# Input:a set of genes of interest, could be from a system or other parts
# Output:plot showing trend of gene expression similarity vs GO semantic similarity
# Dependencies:

# Loading necessary data
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "GEO-QN-profile-2019-03-24.rds")
myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame


corGenes <- function(A, B, prf) {
  # Calculate pearson correlation between gene expression
  # profiles A and B in prf identified by the gene symbol.
  # A and B can be either gene symbol or index.
  
  r <- cor(prf[A, ], prf[B, ], use = "pairwise.complete.obs")
  return(r)
}

# used for sample codes and tests
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

if (! requireNamespace(seqinr, quietly=TRUE)) {
  install.packages("seqinr")
}
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("GOSemSim", version = "3.8")
BiocManager::install("org.Hs.eg.db", version = "3.8")
# packages for calculating GO similarity
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="MF")
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="MF", computeIC=FALSE)
library("biomaRt")
# library(help = seqinr)       # basic information
# browseVignettes("seqinr")    # available vignettes
# data(package = "seqinr")     # available datasets
#================================

  myFunction <- function(A) {
    # Purpose:
    #     Collecting gene correlations and GO semantic similarity and build relationship model
    #     for each pair of genes to prove the hypothesis
    # Parameters:
    #     A : a set of genes of interets
    # Value:
    #     result: summary of the model withe scatterplot.
    
    # code ...
    mydata <- data.frame(GeneA = A[1],GeneB = A[1],CorGene = 1,scoGO = 1)
    len <- length(A)
    for (i in seq_len(len-1)){
      for (j in (i+1) : len){
        if (i != j){
         a <- A[i]
         b <- A[j]
         cor <- corGenes(a,b,myQNXP)
         scoreGO <- geneSim(a, b, semData=hsGO2, measure="Wang", combine="BMA")$geneSim
         new <- data.frame(GeneA = a,GeneB = b,CorGene = cor,scoGO = scoreGO)
         mydata <- rbind(mydata,new)
        }
      }
    }
      mydata <- mydata[-1,]
      x <- mydata$CorGene
      y <- mydata$scoGO
      plot(x, y, main="Gene similarity VS GO semantic similarity", 
           xlab="Gene correlation", ylab="GO score", pch=19)
      model <- lm(y~x)
      result <- summary(model)
      abline(lm(y~x), col="red")
    
    return(result)
  }
 
  
geneSet <- c("AMBRA1","ATG14","ATP2A1","ATP2A2","ATP2A3")
myFunction(geneSet)
