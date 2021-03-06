
# pathwayEnrich.R
#'
#' \code{PathwayEnrichment} a function that perform a pathway enrichment analysis
#' Return a barplot that show the enriched pathway and its p value(fisher test) that is smaller than
#' the ajusted p value. The null hypothesis for each fisher test: the pathway
#' is not enriched in the given gene set.
#'
#' @param geneSet A list of gene symbol in the system that you want to perform pathway
#'  enrichment analysis.
#'
#' @param reactInfo a dataframe map the reactome pathway id to the reference gene symbol
#'
#' @return a dataframe that show the enriched pathway, adjusted pa value using
#' Bonferroni and Benjamini-Hochberg methods,and its p value(fisher test).
#'
#' @import ggplot2
#' @importFrom stats fisher.test
#' @importFrom stats p.adjust
#' @author {Yufei Yang} (aut)
#' @examples
#' \dontrun{
#' # Picking sample gene list to find the pathway enrichment and p value
#' # Call the symComp helper function to generate list
#' sampleGeneList <- c("CREBBP","NR1H3","RORA","SREBF1","SMARCD3")
#' reactInfo <- fetchData("ReactomeSym")
#' result <- PathwayEnrichment(sampleGeneList,reactInfo)
#' }
#' @export
PathwayEnrichment <- function(geneSet, reactInfo){

  #fetch HGNC reference genes
  pathway_df <- data.frame(pathway=c(),pval=c())
  unique_pw <- unique(reactInfo$REACTOME_ID)
  #get data information from the helper function
  fisher_matrix <- matrix(nrow = 2, ncol = 2)

  pwList <- c()
  #get all pathways in the given gene set
  for(gene in geneSet){
    pws <- reactInfo$REACTOME_ID[reactInfo$HGNC==gene]
    for(pw in pws){
      if(!is.na(pw)){
        pwList <- c(pwList,pw)
      }

    }
  }

  uniquePW <- unique(pwList)
  #occurrence list of pathway for given gene set initialized
  occurenceGene <- list()
  #occurrence list of pathway for all human gene initialized
  occurence <- list()
  for(pw in uniquePW){
    occurenceGene[[pw]] <- 0
    occurence[[pw]] <- 0
  }
  #calculate value for the occurrence list of pathway for the given gene set
  for(gene in geneSet){
    pws <- reactInfo$REACTOME_ID[reactInfo$HGNC==gene]
    for(pw in pws){
      if(!is.na(pw)){
        occurenceGene[[pw]] <- occurenceGene[[pw]]+1
      }

    }
  }
  enrichment <- data.frame(pathway=c(),pval=c(),adjustedPvalBon=c())
  #update fisher test matrix
  #calculate value for the occurrence list of pathway for all human gene and given gene set
  for( pathway in pwList){
    #get occurrence from the reactome data
    occur <- sum(reactInfo$REACTOME_ID==pathway)
    occurence[[pathway]] <- occur
    geneList <- reactInfo$HGNC[reactInfo$REACTOME_ID==pathway]
    #fisher matrix column for the all human gene
    #how many time a pw occur in the ref gene
    fisher_matrix[1,2] <- occur
    #how many time this pw are not occur in the ref gene
    fisher_matrix[2,2] <- length(reactInfo$REACTOME_ID) -occur
    occurSet <- occurenceGene[[pathway]]
    #fisher matrix column for the given gene set
    #how many pws are in system
    fisher_matrix[1,1] <- occurSet
    #how many time this pw are not occur in system
    fisher_matrix[2,1] <- length(pwList)-occurSet
    pval <- stats::fisher.test(fisher_matrix,workspace = 200000,alternative = "two.side")$p.value
    temp <- data.frame(pathway=pathway,pval=pval)
    enrichment <- rbind(enrichment,temp)
  }

  #sort by p value get the p value that is larger than 0.05
  #since the null hypothesis is the pathway is enriched
  enrichment <- enrichment[order(enrichment$pval),]
  #get adjusted p value for 2 methods
  len <- length(reactInfo$REACTOME_ID)
  #calculate adjusted p value using BH and bonferroni method
  adjustedPvalBH <- stats::p.adjust(enrichment$pval, method = "BH", n = len)
  adjustedPvalBon <- stats::p.adjust(enrichment$pval, method = "bonferroni", n = len)
  padjust <- data.frame(adjustedPvalBH=adjustedPvalBH,adjustedPvalBon=adjustedPvalBon)
  enrichment <- data.frame(cbind(enrichment,padjust))
  #remove adjusted p value that is equal to 1 (meaning no effect)
  enrichment <- enrichment[enrichment$adjustedPvalBon<1 & enrichment$adjustedPvalBH <1,]
  return(enrichment)

}


# [END]
