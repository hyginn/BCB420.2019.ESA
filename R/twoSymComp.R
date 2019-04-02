# twoSymComp.R

#'
#' \code{twoSymComp} Return a list of domains that the genes in two systems
#'    overlap along with genes that mapped to the same domains
#'
#' @param system1 (character) The name of the system in interest. Is compared with system2
#'     in this function.
#' @param system2 (character) The name of the system in interest. Is compared with system1
#'     in this function.
#' @return (list) A list of overlapping Pfam domains along with genes that mapped to them.
#'     Also, a horizonal bar graph is generated to summarize the count of genes
#'     mapped to each overlapping domains.
#' @author {Tina Lee} (aut)
#' @examples
#' # Picking sample system1 and system2 to get intersect domains of the
#' # two systems
#' # Call the twoSymComp helper function to generate list
#' result <- twoSymComp("PHALY", "SLIGR")
#'
#' @export

twoSymComp <- function(system1, system2) {

  # Load required Pfam domain data
  genesIPR <- fetchData("genesIPR")

  # Get genes in each system
  myDB <- fetchData("SysDB")
  geneSet1 <- SyDBgetSysSymbols(myDB, system1)
  geneSet2 <- SyDBgetSysSymbols(myDB, system2)


  # Map each gene to its correponding domains
  s1 <- genesIPR[unlist(geneSet1)]
  s1 <- s1[!sapply(s1, is.null)]
  s2 <- genesIPR[unlist(geneSet2)]
  s2 <- s2[!sapply(s2, is.null)]

  # find overlapping domains of the two systems
  tmp <- vector()
  for (i in seq_along(s1)) {
    for (j in seq_along(s2)) {
      ints <- unlist(intersect(s1[[i]], s2[[j]]))
      if (length(ints > 0)) {
        tmp <- unique(c(tmp, ints))
      }
    }
  }

  # find corresponding genes that maps to the overlapping domains
  dom <- list()
  for (i in seq_along(tmp)) {
    d <- vector()
    for (j in seq_along(s1)) {
      if (any(s1[[j]] %in% tmp[i]) & (!(names(s1[j]) %in% d))) {
        d <- c(d, names(s1[j]))
      }
    }
    for (k in seq_along(s2)) {
      if (any(s2[[k]] %in% tmp[i]) & (!(names(s2[k]) %in% d))) {
        d <- c(d, names(s2[k]))
      }
    }
    if (length(d) > 0) {
      dom[[i]] <- d
    }
  }
  names(dom) <- tmp

  # Plot horizontal graph of overlapping domains
  res <- sort(unlist(lapply(dom, function(x) length(x))))

  graphics::barplot(res,
                    main = "Overlapping Domains",
                    xlab = "Counts",
                    horiz = TRUE,
                    names.arg = names(res),
                    cex.names = 0.6,
                    col = grDevices::colorRampPalette(c("#FF6655", "#8888A6"),
                                                      bias = 1.3)(length(dom)),
                    las = 1)

  # return a list of overlapping domains along with genes that maps to the
  # domains
  dom <- dom[order(unlist(lapply(dom, function(x) length(x))), decreasing = TRUE)]

  return(dom)
}

# [END]
