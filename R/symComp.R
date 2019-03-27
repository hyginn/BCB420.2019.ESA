# symComp.R
#
# Purpose: Compare domains of two similar systems
# Version: 1.0.0
# Date: March 24, 2019
# Author: Tina Lee (yijen.lee@mail.utoronto.ca)
# License: MIT
#
# Input: Two set of genes from similar systems.
# Output:
# Dependencies: BCB420.2019.Pfam

# ======= Functions =================
# stub function taken from Dr. Steipe's BCB420.2019.ESA package, but changed to
# system I want to work on.
fetchComponents <- function(sys) {
  # returns a fixed set of symbols.
  # Function stub for development purposes only.
  # Got genes for each antigen processing system from KEGG.
  if (sys == "ANPRC") {
    s <- c("IFNG", "TNF", "TNFA", "PSME1", "HSP90A", "HSPA5", "TAP2", " TAP1",
           "TAPBP", "PDIA3", "CANX", "CALR", "B2M", "CD8A", "TRAV", "KLRD1",
           "HLA-A", "HLA-B", "HLA-C", "CIITA", "RFX5", "CREB1", "NFYA", "CD74")
  } else if (sys == "ANPRC2"){
    s <- c("IFI30", "LGMN", "CTSB", "CD74", "HLA-DMA", "CIITA", "RFX5", "CREB1",
           "NFYA", "CD4", "TRAV")
  } else {
    s <- ""
  }
  return(s)
}

# symComp.R
#'
#' \code{symComp} Return a list of domains that the genes in two systems
#'    overlap along with genes that mapped to the same domains
#'
#' @param system1 The name of the system in interest. Is compared with system2
#'     in this function.
#' @param system2 The name of the system in interest. Is compared with system1
#'     in this function.
#' @return A list of overlapping Pfam domains along with genes that mapped to them.
#'     Also, a horizonal bar graph is generated to summarize the count of genes
#'     mapped to each overlapping domains.
#' @author {Tina Lee} (aut)
#' @examples
#' # Picking sample system1 and system2 to get intersect domains of the
#' two systems
#' # Call the symComp helper function to generate list
#' result <- symComp("ANPRC", "ANPRC2")
#'
#' @export

symComp <- function(system1, system2) {

  # Load required Pfam domain data
  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "genesIPR.rds")
  genesIPR <- readRDS(url(myURL))

  # Get genes in each system using fetchComponents()
  geneSet1 <- fetchComponents(system1)
  geneSet2 <- fetchComponents(system2)

  # Map each gene to its correponding domains
  s1 <- genesIPR[geneSet1]
  s1 <- s1[!sapply(s1, is.null)]
  s2 <- genesIPR[geneSet2]
  s2 <- s2[!sapply(s2, is.null)]

  tmp <- vector()
  for (i in seq_along(s1)) {
    for (j in seq_along(s2)) {
      ints <- unlist(intersect(s1[[i]], s2[[j]]))
      if (length(ints > 0)) {
        tmp <- unique(c(tmp, ints))
      }
    }
  }

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

  # Plot horixontal graph of results
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

  # return a list of overlapped domains along with genes that maps to the
  # domains
  dom <- dom[order(unlist(lapply(dom, function(x) length(x))), decreasing = TRUE)]

  return(dom)
}

# [END]
