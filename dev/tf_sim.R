# tf_sim.R

#' Transcription Factor Similarity
#'
#' \code{tf_sim} Calculate pairwise upstream transcription factor binding site similarity.
#'
#' Calculate the similarity of upstream transcription factors of 2 genes.
#' Similarity is calculated as the fraction of TF's of the gene with the most TF's that are shared with the other gene. Thus the maximum possible similarity is 1, achieved when the two TF sets are identical and the minimum, 0, is achieved when no TFs are shared. The similarity metric is independent of numner of TF binding sites.
#' @section <title>: Additional explanation.
#'
#' @param gene1 String, HGNC symbol for the first gene.
#' @param gene2 String, HGNC symbol for the second gene.
#' @param GTRD A list of lists of transcription factors than have binding sites upstream of genes. Loaded from "http://steipe.biochemistry.utoronto.ca/abc/assets/geneList-2019-03-13.RData"
#' @return The similarity score of gene1 and gene2: number of shared upstream transcription factors.
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' # Calculate the similarity of 2 related genes "BRCA1" and "BRCA2"
#' tf_sim("BRCA1", "BRCA2")
#' [1] 0.3809524
#'
#' # calculate similarity of 2 genes that are not known to be related
#' tf_sim("BRCA1", "GDI1")
#' [1] 0.04761905
#'
#' @export

tf_sim <- function(gene1, gene2, GTRD=geneList) {
  tfs1 <- GTRD[[gene1]]
  len1 <- length(tfs1)
  tfs2 <- GTRD[[gene2]]
  len2 <- length(tfs2)
  int <- intersect(tfs1, tfs2)
  similarity <- length(int)/max(c(len1, len2))
  return(similarity)
}

# [END]
