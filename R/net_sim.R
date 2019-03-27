# net_sim.R

#' Network Similarity.
#'
#' \code{net_sim} Calculate the pairwise similarity of genes in an interaction network.
#'
#' Jaccard similarity of 2 nodes in the gene network is calculated. Jaccard similarity is the number of common neighbors divided by the number of vertices that are neighbors of at least one of the 2 nodes.
#' Note that this functions requires the igraph package to previously be loaded.
#' @section <title>: additional explanation
#'
#' @param gene1 Character, HGNC symbol of the first gene.
#' @param gene2 Character, HGNC sysmbol of the second gene
#' @param graph igraph object representing the network
#' @return Numeric, similarity score: Jaccard similarity
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' # Calculate the network similarity of 2 related genes:
#' net_sim("BRCA1", "BRCA2")
#' [1] 0.3481481
#'
#' # Calculate network similarity of 2 genes that are not known to be related:
#' net_sim("ROBO1", "BRCA1")
#' [1] 0.01592357
#'
#' @export

net_sim <- function(gene1, gene2, graph = STRINGgraph) {
  # if the 2 genes are the same and they are unconnected, Jaccard similarity will be 0 but we want it to be 1
  if (gene1==gene2) {
    return(1)
  } else { # compute Jaccard similarity normally
    sim <- similarity(graph, vids = c(gene1, gene2) , method = "jaccard")
    sim <- sim[1, 2]
    return(sim)
  }
}

# [END]
