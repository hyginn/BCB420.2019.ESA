# hypothesize.R
#' @title hypothesize
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' \code{hypothesize} generates an annotated hypothesis graph whose nodes are
#' components of \code{mySys} and edges are GGI interpretations in either
#' \code{EMAP} or \code{GMAP}
#'
#' @param network A dataframe of physically-interacting system components with
#' their genetic interactions and
#' @param ppi_ggi An optional dataframe to specify subset of \code{mySys} for
#' which both PPI and GGI data is available
#'
#' @importFrom dplyr inner_join
#' @import visNetwork
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

#' @title visualizeInteractions
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' A helper function of \code{hypothesize} that plots the hypothesis graph
#'
#' @inheritParams hypothesize
#' @param gmap The dataframe obtained from \code{makeGMAP}.
#' @param emap The dataframe obtained from \code{makeEMAP}.
#'
visualizeInteractions <- function(network, emap, ppi_ggi, gmap) {
  # get all genes to be plotted
  genesInNetwork <- unique(c(as.character(network$A), as.character(network$B)))
  allgenes <- as.factor(genesInNetwork)

  if (is.null(ppi_ggi)) # does not highlight stringent subnetworks
  {
    # prepare nodes and edges information for plotting
    sel <- match(as.character(network$type),
                 as.character(emap$geneticInt))
    edgeLabels <- as.character(emap$effect[sel])
    nodes <- data.frame(id = allgenes, allgenes,
                        label = genesInNetwork,
                        shape = 'circle')
    edges <- data.frame(
        from = network$A,
        to = network$B,
        label = edgeLabels,
        arrows = "to"
      )
  } else {
    # get all genes in stringent subnetwork
    sel <- dplyr::inner_join(network, ppi_ggi)
    genesInPGI <- unique(c(ppi_ggi$A, ppi_ggi$B))
    allgenes2 <- as.factor(genesInPGI)

    # prepare nodes and edges information for plotting
    sel2 <- allgenes %in% allgenes2 & allgenes %in% allgenes2
    groups <- ifelse(sel2 , "ppi-ggi", "ggi")

    nodes <- data.frame(id = allgenes, allgenes, label = as.character(allgenes),
                        group = groups,
                        shape = 'circle')
    sel2 <- network$A %in% ppi_ggi$A & network$B %in% ppi_ggi$B
    gMappings <- match(as.character(network$type),
                      as.character(gmap$geneticInt))
    eMappings <- match(as.character(ppi_ggi$type),
                       as.character(emap$geneticInt))
    edgeLabels <- ifelse(sel2,
                         as.character(gmap$effect[gMappings]),
                        as.character(emap$effect[eMappings]) )

    edges <- data.frame(from = as.factor(network$A),
                        to = as.factor(network$B),
                        label = edgeLabels,
                        arrows = "to")
    }
  visNetwork(nodes, edges) %>%
    visGroups(groupname = "ppi-ggi", color = "salmon") %>%
    visEdges(color = "darkgray")
}

# [END]
