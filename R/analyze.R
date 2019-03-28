#' Function corrGraphs().
# stub function taken from Dr. Steipe's BCB420.2019.ESA package
# Taken from https://github.com/hyginn/
exProf <- function(sym, hgnc, ncol = 20) {
  # returns a set of numbers as a virtual expression profile, for
  # development purposes only.
  set.seed(which(hgnc$sym == sym))
  p <- as.vector(scale(runif(ncol)))
  set.seed(NULL)
  return(p)
}
# stub function taken from Dr. Steipe's BCB420.2019.ESA package
# Taken from https://github.com/hyginn/
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
  } else {
    s <- ""
  }
  return(s)
}


#'
#' \code{corrGraphs()} Calculates Semantic similarity and Correlation
#' between expression profiles of pairwise
#' genes and returns a list of two heat maps representing the two correlation
#' levels and one plot to charcaterize the orthogonal relationship
#' between genes in the system.
#'
#'
#' @param bio.sys Biological system symbol.
#' @return A list of 3 diagrams:
#' l["co-exp"] - Co-expression correlations heatmap;
#' l["semsim"] - Semantic similarity correlations heatmap and
#' l["func"] - Functional similarities
#'
#' @author \href{https://orcid.org/0000-0002-8778-6442}{Denitsa Vasileva} (aut)
#'
#' @examples
#' corrGraphs("PHALY")
#'
#' @export
corrGraphs <- function(bio.sys) {
  # following segment is adapted from Dr. Steipe's BCB420.2019.ESA package
  myURL <- paste0("https://github.com/hyginn/",
                  "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")

  load(url(myURL))  # loads HGNC data frame

  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "geneList-2019-03-13.RData")

  load(url(myURL))

  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "STRINGedges-2019-03-14.RData")
  load(url(myURL))  # loads STRING edges object

  myURL <- paste0("https://github.com/hyginn/",
                  "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))  # loads HGNC data frame

  myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                  "GEO-QN-profile-2019-03-24.rds")
  myQNXP <- readRDS(url(myURL))
  #end of adapted segment

  cat(sprintf("Fetching components for %s...\n",bio.sys))
  system.components <- fetchComponents(bio.sys)
  if (system.components == "")
    stop("Unknown system passed as parameter.")
  num_components <- length(system.components)
  cat(sprintf("%d components found for %s.\n",num_components, bio.sys))

  cat(sprintf("Loading GO - genome wide annotation for human.\n"))
  hsGO <- GOSemSim::godata('org.Hs.eg.db', keytype = "SYMBOL", ont = "MF")
  cat(sprintf("%d rows loaded.\n",nrow(hsGO)))

  # calculate correlations
  gene1 <- vector(mode = "character", length = 0)
  gene2 <- vector(mode = "character", length = 0)
  correlation <- vector(mode = "numeric", length = 0)
  go_correlation <- vector(mode = "numeric", length = 0)
  corrM <- matrix(nrow = num_components, ncol = num_components)
  goM <- matrix(nrow = num_components, ncol = num_components)
  colnames(corrM) <- system.components
  rownames(corrM) <- system.components
  colnames(goM) <- system.components
  rownames(goM) <- system.components

  cat(sprintf("Calculating correlations..."))
  for (component1 in 1:num_components) {
    for (component2 in component1:num_components) {
      prf1 <- myQNXP[system.components[component1],]
      prf1[is.na(prf1)] <- 0
      prf2 <- myQNXP[system.components[component2],]
      prf2[is.na(prf2)] <- 0
      myQNXP[system.components[component1],]
      cr1 <- cor(prf1, prf2, use = "pairwise.complete.obs")
      cr2 <- GOSemSim::geneSim(system.components[component1],
                               system.components[component2],
                               semData = hsGO,
                               measure = "Wang",
                               combine = "BMA")[[1]]
      gene1 <- c(gene1, system.components[component1])
      gene2 <- c(gene2, system.components[component2])
      correlation <- c(correlation, cr1)
      go_correlation <- c(go_correlation, cr2)
      corrM[component1, component2] <- cr1
      corrM[component2, component1] <- cr1
      goM[component1, component2]   <- cr2
      goM[component2, component1]   <- cr2
    }
    cat(sprintf("."))
  }
  cat(sprintf("\n"))

  g.corr <- corrplot::corrplot(corrM,
                     type = "upper",
                     order = "original",
                     tl.cex = 0.5,
                     diag = TRUE,
                     title = "Functional Correlation",
                     col = RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")
  )
  g.semsim <- corrplot::corrplot(goM,
                     type = "upper",
                     order = "original",
                     tl.cex = 0.5,
                     diag = TRUE,
                     na.label = "N",
                     na.label.col = "white",
                     title = "Semantic similarity correlation",
                     col = RColorBrewer::brewer.pal(n = 8, name = "RdYlBu")
  )

  df <-  data.frame(gene1, gene2, correlation, go_correlation)
  df <- na.omit(df)
  gpl <- ggplot2::ggplot(df, ggplot2::aes(gene1, gene2 )) +
    ggplot2::geom_tile(ggplot2::aes(fill = correlation), , color = "white") +
    ggplot2::scale_fill_gradient(low = "red", high = "green") +
    #ggplot2::ylab("List of genes") +
    #ggplot2::xlab("List of genes") +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 10),
          legend.text = ggplot2::element_text(size = 12),
          plot.title = ggplot2::element_text(size=16),
          axis.title = ggplot2::element_text(size=14,face="bold"),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(fill = "Correlation")

  cat(sprintf("Functional and co-expression correlation diagram\n"))
  g.func <- ggplot2::ggplot(data = df,
                ggplot2::aes(x = go_correlation, y = correlation)) +
                ggplot2::geom_point(shape = 22, size = 2.5, colour = "Red") +
                ggplot2::ggtitle("Functional and co-expression correlation") +
                ggplot2::scale_x_continuous(name = "Semantic similarity") +
                ggplot2::scale_y_continuous(name = "Co-expression correlation")
  ## + geom_text(aes(label = title)) # should you need labels

  return(list("co-exp" = g.corr,
              "semsim" = g.semsim,
              "func"   = g.func))
}
