# fetchComponents() adapted from BCB420.2019.ESA by Dr. Boris Steipe
# https://github.com/hyginn/BCB420.2019.ESA

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
  }
  # Add high voltage calcium channel system curated by Edward Ho
  else if (sys == "VGCR") {
    s <- c("AKAP7", "PRKACA", "PRKACB", "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B",
           "ADCY9", "ADRB2", "GNAS", "CALM1", "CALM2", "CALM3", "CACNA1C", "CACNA2D1",
           "CACNA2D3", "CACNB1", "CACNB3", "CACNB4", "CACNG1", "CACNG6")
  }

  else {
    s <- ""
  }
  return(s)
}


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
  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "STRINGedges-2019-03-14.RData")
  load(url(myURL))

  # Create graph from STRING edges
  geneGraph <- igraph::graph_from_edgelist(matrix(c(STRINGedges$a, STRINGedges$b), ncol = 2, byrow = FALSE), directed = FALSE)

  # Get genes that are system components
  geneSet <- fetchComponents(sys)

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
