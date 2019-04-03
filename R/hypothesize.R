# hypothesize.R
#' @title hypothesize
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' \code{hypothesize} generates an annotated hypothesis graph whose nodes are
#' components of \code{mySys} and edges are GGI interpretations in either
#' \code{EMAP} or \code{GMAP}
#'
#' @param network (dataframe) 3-column dataframe of physically-interacting
#' system components with their genetic interactions. The first two columns
#' denote the interacting pair; the third is the type of genetic interaction.
#' @param ppi_ggi (dataframe) An optional 3-column dataframe to specify subset
#' of \code{mySys} for which both PPI and GGI data is available. The first two
#' columns denote the interacting pair; the third is the type of genetic
#' interaction.
#' @return (NULL) The function plots the graph of system components.
#'
#' @importFrom dplyr inner_join
#' @import visNetwork
#' @examples
#' # Plot the graph of SLIGR components under stringent and relaxed conditions
#' mySys <- getSysInteractions("SLIGR", criterion = "stringent")
#' mySys2 <- getSysInteractions("SLIGR", criterion = "relaxed")
#' hypothesize(mySys) # draws hypothesis graph(s)
#' hypothesize(mySys2, mySys) # draws hypothesis graph(s)
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
#' @param gmap (dataframe) An 11-by-2 dataframe mapping BioGrid GGI tag to its
#' interpretation. The first column contains the official BioGRID GGI tags; the
#' second contains interpreted relationships.
#' @param emap T(dataframe) An 11-by-3 dataframe mapping BioGrid GGI tag to its
#'  interpretation assuming that the system's components also interact
#'  physically. The first column contains the official BioGRID GGI tags; the
#'  second contains interpreted relationships under the assumption; the third
#'   contains notes on interpretation.
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
                        as.character(emap$effect[eMappings]),
                        as.character(gmap$effect[gMappings]))

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
