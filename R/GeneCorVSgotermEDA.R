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


# ====  PACKAGES  ==============================================================
# Load all required packages.
#
# Use non-standard libraries with  package::function() idiom if possible.


#BiocManager::install("GOSemSim", version = "3.8")
#BiocManager::install("org.Hs.eg.db", version = "3.8")
#install.packages("ggcorrplot")
# packages for calculating GO similarity




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

#geneSet <- c("AMBRA1","ATG14","ATP2A1","ATP2A2","ATP2A3")
#myPlot(geneSet)


# [END]
