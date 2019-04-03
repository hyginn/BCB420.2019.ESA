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
# Dependencies: devtools,ggplot,ggcorrplot, GOSemSim,org.Hs.eg.db
#

# ====  FUNCTIONS  =============================================================


#' \code{installHumanGenomeAnnotation} install Human genome Annotation when necessary.
#'

#' @return (NULL)
#' #' @author {Yuhan Zhang} (aut)
#' @examples
#' \dontrun{
#' installHumanGenomeAnnotation()
#' }
#' @export
installHumanGenomeAnnotation <- function() {
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("org.Hs.eg.db", version = "3.8")
return(invisible(NULL))
}

# function stub taken from Dr. Steipe's BCB420.2019.ESA package (https://github.com/hyginn/BCB420.2019.ESA)
corGenes <- function(A, B, prf) {
  # Calculate pearson correlation between gene expression
  # profiles A and B in prf identified by the gene symbol.
  # A and B can be either gene symbol or index.

  r <- cor(prf[A, ], prf[B, ], use = "pairwise.complete.obs")
  return(r)
}



#'
#' \code{myScatterPlot()}
#' Investigates relationship of expression correlation and semantic similarity of input genes
#' and produces plots and summary of model info
#'
#' @param geneSet A set of genes of interest to investigate
#' @return (list) list of graphs containing A heatmap of co-expression correlation of each pair of genes in the geneSet,
#'                  A heatmap of semantic similarity of each pair of genes,and
#'                  A scatterplot of correlation vs Go similarity with summary of linear model built by correlation vs Go similarity
#'
#' [GOSemSim package](http://bioconductor.org/packages/release/bioc/html/GOSemSim.html)
#' [org.Hs.eg.db data package](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

#'
#' @author {Yuhan Zhang} (aut)
#'
#' @examples
#' geneSet <- c("AMBRA1","ATG14","ATP2A1","ATP2A2","ATP2A3")
#' myScatterPlot(geneSet)
#'
#' @export
  myScatterPlot <- function(geneSet) {
    # Purpose:
    #     Collecting gene correlations and GO semantic similarity and build relationship model
    #     for each pair of genes to prove the hypothesis
    # Parameters:
    #     A : a set of genes of interets
    # Value:
    #     result: summary of the model withe scatterplot.

    # code ...
    # Loading necessary data

    myQNXP <- fetchData("GEOprofiles")  # loads quantile-normalized expression data
    HGNC <- fetchData("HGNCreference")
    installHumanGenomeAnnotation()
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

    geneDataFrame <- data.frame(GeneA, GeneB, CorGene,scoGO)
    corM <- matrix(nrow = len, ncol = len)
    goM <- matrix(nrow = len,ncol = len)

    colnames(corM) <- unique(geneDataFrame$GeneB)
    rownames(corM) <- unique(geneDataFrame$GeneA)
    colnames(goM) <- unique(geneDataFrame$GeneB)
    rownames(goM) <- unique(geneDataFrame$GeneA)

    for (i in 1:nrow(geneDataFrame)) {
      corM[
        geneDataFrame$GeneA[i],
        geneDataFrame$GeneB[i]
        ] <- CorGene[i]
      corM[
        geneDataFrame$GeneB[i],
        geneDataFrame$GeneA[i]
        ] <- CorGene[i]
      goM[
        geneDataFrame$GeneA[i],
        geneDataFrame$GeneB[i]
        ] <- scoGO[i]
      goM[
        geneDataFrame$GeneB[i],
        geneDataFrame$GeneA[i]
        ] <- scoGO[i]
    }
    corM[is.na(corM)] <- 0
    heatmap1 <- ggcorrplot::ggcorrplot(corM, hc.order = TRUE, type = "lower",
                           outline.col = "white", lab = TRUE,title = "Expression Correlation")
    heatmap2 <- ggcorrplot::ggcorrplot(goM, hc.order = TRUE, type = "lower",
                           outline.col = "white", lab = TRUE, title = "GO Term Similarity")
    #print(heatmap1)
    #print(heatmap2)
      pairFrame <- geneDataFrame[which(geneDataFrame$GeneB != geneDataFrame$GeneA),]
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
     #print(lmplot)
    print(summary(model))
    result <- list(heatmap1, heatmap2,lmplot)

    return(result)
  }

#geneSet <- c("AMBRA1","ATG14","ATP2A1","ATP2A2","ATP2A3")
#myScatterPlot(geneSet)


# [END]
