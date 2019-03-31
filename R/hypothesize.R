# hypothesize.R

#' Generates an annotated hypothesis graph whose nodes are components of \code{mySys} and edges are GGI interpretations in either \code{EMAP} or \code{GMAP}
#'
#' @param network A dataframe of physically-interacting system components with their genetic interactions and
#' @param ppi_ggi An optional dataframe to specify subset of \code{mySys} for which both PPI and GGI data is available
#' @import visNetwork
#' @import xlsx
#' @import readxl
#' @importFrom dplyr filter
#' @import biomaRt
#' @import ggplot2
#' @import biogridr
#' @import utils
#' @include SyDButils.R fetchComponents.R fetchData.R
#' @importFrom stats complete.cases
#' @examples
#' myKey <- getKey("Nada Elnour", "nada.elnour@@mail.utoronto.ca", "SILGRESA")
#' mySys <- getSysInteractions("SLIGR", criterion = "stringent", key = myKey)
#' mySys2 <- getSysInteractions("SLIGR", criterion = "relaxed", key = myKey)
#' hypothesize(mySys)
#' hypothesize(mySys2, mySys)
#'
#' @export
hypothesize <-
  function(network,
           ppi_ggi = NULL) {
    EMAP <- makeEMAP()
    GMAP <- makeGMAP()
    visualizeInteractions(network, EMAP, ppi_ggi, GMAP)
  }

#' A helper function of \code{hypothesize} that plots the hypothesis graph
#'
#' @inheritParams hypothesize
#' @param gmap The dataframe obtained from \code{makeGMAP}.
#' @param emap The dataframe obtained from \code{makeEMAP}.
visualizeInteractions <- function(network, emap, ppi_ggi, gmap) {
  if (is.null(ppi_ggi))
  {
    network$gene1 <- as.factor(network$gene1)
    network$gene2 <- as.factor(network$gene2)
    allgenes <-
      as.factor(unique(c(
        as.character(network$gene1),
        as.character(network$gene2)
      )))
    nodes <-
      data.frame(id = allgenes,
                 allgenes,
                 label = unique(c(
                   as.character(network$gene1),
                   as.character(network$gene2)
                 )),
                 shape = 'circle')
    edges <-
      data.frame(
        from = network$gene1,
        to = network$gene2,
        label = emap$effect[match(network$interactionType, emap$geneticInt)],
        arrows = "to"
      )
  } else {
    network$gene1 <- as.factor(network$gene1)
    network$gene2 <- as.factor(network$gene2)

    ppi_ggi$gene1 <- as.factor(ppi_ggi$gene1)
    ppi_ggi$gene2 <- as.factor(ppi_ggi$gene2)

    allgenes <-
      as.factor(unique(c(
        as.character(network$gene1),
        as.character(network$gene2)
      )))

    allgenes2 <-
      as.factor(unique(c(
        as.character(ppi_ggi$gene1),
        as.character(ppi_ggi$gene2)
      )))

    nodes <-
      data.frame(
        id = allgenes,
        allgenes,
        label = unique(c(
          as.character(network$gene1),
          as.character(network$gene2)
        )),
        group = ifelse(allgenes %in% allgenes2 , "ppi-ggi", "ggi"),
        shape = 'circle'
      )

    edgeLabels <- c()

    for (i in 1:(length(network$interactionType))) {
      sel <-
        (
          as.character(network$gene1[i]) %in% as.character(ppi_ggi$gene1) &
            as.character(network$gene2[i]) %in% as.character(ppi_ggi$gene2)
        )
      if (sel) {
        idx <-
          which(
            as.character(ppi_ggi$gene1) == as.character(network$gene1[i]) &
              as.character(ppi_ggi$gene2) == as.character(network$gene2[i])
          )
        edgeLabels <-
          c(edgeLabels, as.character(emap$effect[match(ppi_ggi$interactionType[idx], emap$geneticInt)]))
      } else {
        edgeLabels <-
          c(edgeLabels, as.character(gmap$effect[match(network$interactionType[i], gmap$geneticInt)]))
      }
    }

    edges <-
      data.frame(
        from = network$gene1,
        to = network$gene2,
        label = edgeLabels,
        arrows = "to"
      )

  }

  visNetwork(nodes, edges) %>% visGroups(groupname = "ppi-ggi", color = "salmon") %>% visOptions(selectedBy = "group") %>% visEdges(color = "darkgray")
}

# [END]
