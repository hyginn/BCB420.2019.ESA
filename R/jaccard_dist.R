# jaccard_dist.R

#' Jaccard Distance.
#'
#' \code{jaccard_dist} calculates the pairwise distance of genes in an interaction network.
#'
#' Jaccard similarity of 2 nodes in the gene network is calculated. Jaccard similarity
#' is the number of common neighbors divided by the number of vertices that are
#' neighbors of at least one of the 2 nodes. Distance is then obtained by taking
#' 1 - similarity.
#'
#' @param gene1 Character, HGNC symbol of the first gene.
#' @param gene2 Character, HGNC sysmbol of the second gene
#' @param graph igraph object representing the network
#' @return Numeric, distance score
#'
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#' @import igraph
#'
#' @examples
#'
#' STRING <- fetchData("STRINGedges0.8")
#' # convert the STRINGedges object into an igraph object
#' if (requireNamespace("igraph")) {
#'     STRINGgraph <- igraph::graph_from_edgelist(as.matrix(STRING[,1:2]))
#' } else {
#'     cat("igraph is required for this example")
#' }
#'
#' # calculate the Jaccard network distance between "BRCA1" and "BRCA2"
#' jaccard_dist("BRCA1", "BRCA2", STRINGgraph)
#'
#' @export

jaccard_dist <- function(gene1, gene2, graph) {
  # if the 2 genes are the same and they are unconnected, Jaccard similarity will
  # be 0 but we want it to be 1
  if (gene1 == gene2) {
    sim <- 1
  } else { # compute Jaccard similarity normally
    sim <- similarity(graph, vids = c(gene1, gene2) , method = "jaccard")
    sim <- sim[1, 2]
  }
  distance <- 1 - sim
  return(distance)
}

# [END]
