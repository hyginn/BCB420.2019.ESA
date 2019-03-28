# findSTRINGKnodes.R

#' \code{findSTRINGKnodes} output a the K node values for vertices in the STRING interaction graph given the system of genes 'sys'.
#'
#' \code{findSTRINGKnodes} outputs a the K node values (as defined by Cornish & Markowetz, 2014),
#' a measure of the vertex association to the high-weighted vertices defined by the input
#' curated gene system 'sys'. Also plots the knode values for the top percentile of knode values
#' given by 'per'
#'
#' @param sys (character)       5-character system code of genes
#' @param per (numeric)         Percentile of results to be plotted. Default is 0.05

#' @return (numeric) Knode values
#' @examples
#' findSTRINGKnodes("PHALY", 0.05)

#' @export
findSTRINGKnodes <- function(sys, per=0.05) {

  # Load high-confidence STRING edges
  STRINGedges <- fetchData("STRINGedges0.9")

  # Create graph from STRING edges
  geneGraph <- igraph::graph_from_edgelist(matrix(c(STRINGedges$a, STRINGedges$b), ncol = 2, byrow = FALSE), directed = FALSE)

  # Get genes that are system components
  myDB <- fetchData("SysDB")
  geneSet <- SyDBgetSysSymbols(myDB, sys)

  # Find IDs of vertices that are part of geneSet
  ids <- match(geneSet, igraph::V(geneGraph)$name)
  ids <- ids[!is.na(ids)]

  # Create a vertex attribute called "annotated" that identifies if it was annotated as part of the
  # input system (sys)
  geneGraph <- igraph::set_vertex_attr(geneGraph, name="annotated", value=0)
  geneGraph <- igraph::set_vertex_attr(geneGraph, name="annotated", index=ids, value=1)

  # Calculate Knodes (rank of vertex association to annotated vertices)
  # https://www.rdocumentation.org/packages/SANTA/versions/2.10.2/topics/Knode
  knodes <- SANTA::Knode(geneGraph, dist.method = c("shortest.paths"), vertex.attr = "annotated", verbose = TRUE)

  plot(knodes[1:floor(length(knodes)*per)])

  return (knodes)
}

# [END]
