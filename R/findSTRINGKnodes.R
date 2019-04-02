# findSTRINGKnodes.R

#' \code{findSTRINGKnodes} output the K node values for vertices in the STRING interaction graph given the system of genes 'sys'.
#'
#' \code{findSTRINGKnodes} outputs a the K node values (as defined by Cornish & Markowetz, 2014),
#' a measure of the vertex association to the high-weighted vertices defined by the input
#' curated gene system 'sys'. Also plots the knode the top \emph{n} knode values
#'
#' @param sys (character)       5-character system code of genes

#' @return (numeric) Knode values
#' @examples
#' # Find the Knodes values for genes in the STRING graph sorted by association to PHALY genes
#' findSTRINGKnodes("PHALY")

#' @export
findSTRINGKnodes <- function(sys) {

  # Load high-confidence STRING edges
  STRINGedges <- fetchData("STRINGedges0.9")

  # Create graph from STRING edges
  geneGraph <- igraph::graph_from_edgelist(matrix(c(STRINGedges$a, STRINGedges$b), ncol = 2, byrow = FALSE), directed = FALSE)

  # Get genes that are system components
  myDB <- fetchData("SysDB")
  geneSet <- SyDBgetSysSymbols(myDB, sys)[[1]]

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

  # Within the first n nodes in the ranked list (knodes), find which genes
  # were not included in geneSet (genes in sys)
  topRankedGenes <- names(knodes[1:length(geneSet)])
  inGeneSetID <- which(topRankedGenes %in% geneSet)
  notInGeneSetID <- which(!topRankedGenes %in% geneSet)

  # Create R dataframe made of the top n nodes in the ranked list for plotting
  numPlotPts <- length(geneSet)
  knodeDF <- data.frame(knodes)
  knodeDF <- knodeDF[1:numPlotPts, , drop=FALSE] # drop parameter to retain row names
  knodeDF$inGeneSet <- topRankedGenes %in% geneSet # Mark the genes that were in the original system sys

  # Plot scatterpoint of above dataframe, and have the genes that were in sys
  # labelled in green, genes that were not in sys labelled in red

  ggplot2::ggplot(knodeDF, ggplot2::aes(x=1:numPlotPts, y=knodes, label = row.names(knodeDF), color = factor(inGeneSet, levels=c(TRUE, FALSE)))) +
    ggrepel::geom_label_repel() +
    ggplot2::geom_point(color = 'black') +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::scale_color_manual(values = c("green", "red"), name = "In gene set") +
    ggplot2::labs(x="Gene Knode Rank", y="Knode Value", title=paste("Top Knode Values Through Association with Genes in ", sys))

  return (knodes)
}

# [END]
