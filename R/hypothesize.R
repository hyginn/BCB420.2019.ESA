# hypothesize.R

#' Generates an annotated hypothesis graph whose nodes are components of \code{mySys} and edges are GGI interpretations in either \code{EMAP} or \code{GMAP}
#'
#' @param network A dataframe of physically-interacting system components with their genetic interactions and
#' @param ppi_ggi An optional dataframe to specify subset of \code{mySys} for which both PPI and GGI data is available
#' @import visNetwork
#'
#' @examples
#' mySys <- getSysInteractions("SLIGR", criterion = "stringent")
#' mySys2 <- getSysInteractions("SLIGR", criterion = "relaxed")
#' hypothesize(mySys)
#' hypothesize(mySys2, mySys)
#'
#' @export
hypothesize <- function(network, ppi_ggi = NULL) {
    # make interpretation maps
    EMAP <- makeEMAP()
    GMAP <- makeGMAP()

    # plot the network
    visualizeInteractions(network, EMAP, ppi_ggi, GMAP)
  }

#' A helper function of \code{hypothesize} that plots the hypothesis graph
#'
#' @inheritParams hypothesize
#' @param gmap The dataframe obtained from \code{makeGMAP}.
#' @param emap The dataframe obtained from \code{makeEMAP}.
visualizeInteractions <- function(network, emap, ppi_ggi, gmap) {
  if (is.null(ppi_ggi)) # if user does not want to highlight stringent subnetworks
  {
    genesInNetwork <- unique(c(network$A, network$B))

    network$A <- as.factor(network$A)
    network$B <- as.factor(network$B)
    allgenes <- as.factor(genesInNetwork)

    # prepare nodes and edges information for plotting
    nodes <- data.frame(id = allgenes, allgenes, label = genesInNetwork, shape = 'circle')
    edges <- data.frame(
        from = network$A,
        to = network$B,
        label = emap$effect[match(network$type, emap$geneticInt)],
        arrows = "to"
      )
  } else {
    genesInNetwork <- unique(c(as.character(network$A), as.character(network$B)))
    allgenes <- as.factor(genesInNetwork)

    genesInPGI <- unique(c(as.character(ppi_ggi$A), as.character(ppi_ggi$B)))
    allgenes2 <- as.factor(genesInPGI)

    network$A <- as.factor(network$A)
    network$B <- as.factor(network$B)

    ppi_ggi$A <- as.factor(ppi_ggi$A)
    ppi_ggi$B <- as.factor(ppi_ggi$B)

    # prepare nodes and edges information for plotting
    nodes <- data.frame(id = allgenes, allgenes, label = genesInNetwork,
                        group = ifelse(allgenes %in% allgenes2 , "ppi-ggi", "ggi"), shape = 'circle')

    edgeLabels <- c()

    for (i in 1:(length(network$type))) {

      sel <- (as.character(network$A[i]) %in% as.character(ppi_ggi$A) &
                as.character(network$B[i]) %in% as.character(ppi_ggi$B))

      if (sel) {
        idx <- which(as.character(ppi_ggi$A) == as.character(network$A[i]) &
                       as.character(ppi_ggi$B) == as.character(network$B[i]))
        interpretation <- emap$effect[match(ppi_ggi$type[idx], emap$geneticInt)]
        edgeLabels <- c(edgeLabels, as.character(interpretation))
      } else {
        interpretation <- gmap$effect[match(network$type[i], gmap$geneticInt)]
        edgeLabels <- c(edgeLabels, as.character(interpretation))
      }
    }
    edges <- data.frame(from = network$A, to = network$B, label = edgeLabels, arrows = "to")
  }

  visNetwork(nodes, edges) %>% visGroups(groupname = "ppi-ggi", color = "salmon") %>% visEdges(color = "darkgray")
}

# [END]
