# subunitcoExpressionAnalysis.R

# Utility functions for complex subunits coexpression analysis
# to discover outliers in a complex. This will give us insight
# about subunits with particular functional roles or
# allow us to detect anomalies of a complex
#   - processmyQNXP()
#   - getComplexmyQNXP()
#   - complexExpHeatmap()
#   - getCoExpProfile()
#   - getpValue()
#   - genePairsShuffle()
#   - pairExpHist()


#' \code{processmyQNXP} Load and process expression-profiles
#'
#' \code{processmyQNXP} The expression profiles were compiled and normalized by
#' using GEO data containing 52 microarray experiments.
#' For details see \code{\link{fetchData}}. The expression profiles were then
#' processed to eliminate genes with no data (row with NAs cross all experiments)
#'
#' @return (matrix) A matrix with 21063 rows representing genes in HGNC symbols
#' and 52 columns representing the 52 microarray experiments. Each cell of the matrix
#' contains the corresponded expression values
#'
#' @seealso \code{\link{fetchData}}
#'
#' @examples
#' \dontrun{
#' # fetch and process expression profile data
#' myQNXPProcesed <- processmyQNXP()
#' }
#'
processmyQNXP <- function(){
  # load GEO data
  myQNXP <- fetchData("GEOprofiles")
  # delete rows with no data across all experiments
  myQNXPProcesed <- myQNXP[rowSums(is.na(myQNXP)) != ncol(myQNXP), ]
  return(myQNXPProcesed)
}

#' \code{getComplexmyQNXP} Extract expression profile for the complex subunits
#'
#' \code{getComplexmyQNXP} The processed expression profile from GEO is extracted
#'  for complex subunits in both CORUM database and GEO. For details of CORUM data
#'  fetch see \code{\link{fetchData}}.
#'
#' @param complex (char)  name of the complex
#'
#' @return (matrix) A matrix with rows representing complex subunits in HGNC symbols
#' and 52 columns representing the 52 microarray experiments. Each cell of the matrix
#' contains the corresponded expression values
#'
#' @seealso \code{\link{fetchData}}
#'
#' @examples
#' \dontrun{
#' complex <- "MDC1-MRN-ATM-FANCD2 complex"
#' # get the expression file of the comples subunits
#' myComplexQNXP <- getComplexmyQNXP(complex)
#' }
#'

getComplexmyQNXP <- function(complex){
  # load CORUM data
  CORUM <- fetchData("complexSymbol")
  # get the subunits of the comples in HGNC symbols
  complexGene <- CORUM[which(names(CORUM) == complex)][[1]]
  cat(paste("The subunits of", complex, ":\n"))
  print(complexGene)
  # extract information for genes in both the CORUM and the GEO
  myQNXP <- processmyQNXP()
  complexGene <- intersect(complexGene, row.names(myQNXP))
  myComplexQNXP <- myQNXP[complexGene, ]
  return(myComplexQNXP)
}


#' \code{complexExpHeatmap} Heatmap of the complex subunits expression profile
#'
#' \code{complexExpHeatmap} Generate a Heatmap of complex subunits expression profile
#' for visualization. Each row for the heatmap displays the gene encoding for each
#' subunit in HGNC symbols, whereas each column displays an experiment. For each grid,
#' the more navy it is, the more negatively expressed the gene is. Vice versa for the
#' brightness of yellow representing positive expression
#'
#' @param complex (char)  name of the complex
#'
#' @return (list) A plot representing expression heatmap of the complex subunits
#'                  across all microarray experiments
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' complex <- "MDC1-MRN-ATM-FANCD2 complex"
#' heatmap <- complexExpHeatmap(complex) # give the expression profile heatmap
#'                                       # for visualization
#' heatmap # display the heatmap
#' }
#' @export
complexExpHeatmap <- function(complex){
  myComplexQNXP <- getComplexmyQNXP(complex)
  # log transform the expression values
  myComplexQNXP <- log(myComplexQNXP + 1, 2)
  myComplexQNXP <- scale(myComplexQNXP)
  # centering the expression values
  myComplexQNXPdf <- as.data.frame(myComplexQNXP)
  myComplexQNXPdf <- tibble::rownames_to_column(myComplexQNXPdf, var = "genes")
  myComplexQNXPdf <- tidyr::gather(myComplexQNXPdf, "experiments", "expValue", -"genes")

  # Generate the expression heatmap
  expHeatmap <- ggplot2::ggplot(myComplexQNXPdf,
                                ggplot2::aes(x = myComplexQNXPdf$experiments, y = myComplexQNXPdf$genes, fill = myComplexQNXPdf$expValue)) +
    ggplot2::geom_tile()+
    ggplot2::scale_fill_gradient2(mid = "white", low="navy", high="yellow", midpoint = 0) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(size = 9, angle = 330,
                                                       hjust = 0, colour = "black")) +
    ggplot2::ggtitle(paste(complex, "Expression Heatmap")) +
    ggplot2::xlab("Conditions") + ggplot2::ylab("Genes") +
    ggplot2::labs(fill = "Expression \n Values") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  return(expHeatmap)
}


#' \code{getCoExpProfile} Generate a coexpression profile
#'
#' \code{getCoExpProfile} Generate a coexpression profile with all possible pairs of
#' genes in the complex and their Pearson Correlation Coefficients represent how
#' coexpressed the gene pairs are. This function also displays a histogram to show the
#' distribution of the coexpression correlation. See \code{\link{EnrichDepletTF}} for
#' expProfile generate (line 149-154)
#'
#' @param complex (char)  name of the complex
#'
#' @return (data frame) a data frame with there columns: the first gene, the second gene
#' and the Pearson Correlation Coefficients for their coexpression
#'
#' @importFrom stats cor
#' @import graphics
#'
#' @seealso \code{\link{EnrichDepletTF}}
#'
#' @examples
#' complex <- "MDC1-MRN-ATM-FANCD2 complex"
#' expProfile <- getCoExpProfile(complex)
#' @export
getCoExpProfile <- function(complex){
  myComplexQNXP <- getComplexmyQNXP(complex)
  complexGene <- row.names(myComplexQNXP)
  # number of possible pairs of genes for the complex
  genePairs <- (length(complexGene) * (length(complexGene) - 1)) / 2
  # initiate a dataframe for coexpression profile data
  expProfile <- as.data.frame(matrix(nrow = genePairs, ncol = 3, data = c("")),
                              stringsAsFactors = FALSE)
  colnames(expProfile) <- c("Gene1", "Gene2", "Correlation")
  # assign the first gene
  expProfile$Gene1 <- unlist(lapply(complexGene,
                                    function (x) rep(x = x,
                                                     (length(complexGene) - which(complexGene == x)))))
  # assign the second genes
  expProfile$Gene2 <- unlist(lapply(complexGene[2 : length(complexGene)],
                                    function (x) complexGene[(which(complexGene == x)) : length(complexGene)]))

  # calculate and assign coexpression correlation for each pair of genes
  for(i in 1:genePairs){
    expProfile$Correlation[i] <- stats::cor(myComplexQNXP[(expProfile$Gene1[i]), ], myComplexQNXP[(expProfile$Gene2[i]), ], use = "pairwise.complete.obs")
  }
  # generate a histogram for distribution visualization of the correlation
  expProfile$Correlation <- as.numeric(expProfile$Correlation)
  graphics::hist(expProfile$Correlation, xlim = c(-1, 1), col = "light blue",
       main = paste(complex, "Expression Profile Correlations"),
       ylab = "Density", xlab = "Coefficient of Correlation [1, -1]")
  return(expProfile)
}

#' \code{getpValue} calculate p-Value
#'
#' \code{getpValue} Calculate the p value of a gene pair correlation on
#' a simulated random gene pair coexpression corelation distribution.
#'
#' @param simulationCorr (double vector)  a vector of correlation grnerate by
#'                                        randomly simulate N gene pairs and
#'                                        calculate correlation.
#'
#' @param myGenePairCorr (double)         gene pair coexpression correlation
#'                                        that we are interested in
#'
#' @return (double) p-Value of how differentially coexpressed the gene pair is
#'                  over the random coexpression correlation
#'
#' @importFrom stats sd pnorm
#'
#' @examples
#' \dontrun{
#' simulationCorr <- c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)
#' myGenePairCorr <- 0.1
#' getpValue(simulationCorr, myGenePairCorr) # 0.7630246
#' }
getpValue <- function(simulationCorr, myGenePairCorr){
  # calculate the standard deviation
  s <- stats::sd(simulationCorr, na.rm = TRUE)
  # calculate the mean
  xbar <- mean(simulationCorr, na.rm = TRUE)
  # calculate the z-scores
  z <- (xbar - myGenePairCorr) / s
  # calculate the p-Value
  pValue <- 2*stats::pnorm(-abs(z))
  return(pValue)
}

#' \code{genePairsShuffle} Generate Coexpression correlation generated by
#' shuffled expression values of the two genes
#'
#' \code{genePairsShuffle} Coexpression correlation created by shuffled
#' expression values of the two genes by 100000 times
#'
#' @param gene1 (char)  HGNC symbol for the first gene of the coexpression
#'                      profile
#' @param gene2 (char)  HGNC symbol for the second gene
#'
#' @return (double vector) A double vector of Coexpression correlation generated
#' by shuffled expression values of the two genes
#'
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' gene1 <- "NBN"
#' gene2 <- "ATM"
#' corrShuffled <- genePairsShuffle(gene1, gene2)
#' corrShuffled
#'}
#'
genePairsShuffle <- function(gene1, gene2){
  N <- 100000
  myComplexQBXP <- getComplexmyQNXP(complex)
  corrShuffled <- c()
  # repeat the for loop for 100000 times
  for (i in 1:N){
    # shuffle expression values of gene 1 and gene2
    gene1Exp <- sample(myComplexQBXP[gene1, ])
    gene2Exp <- sample(myComplexQBXP[gene2, ])
    # calculate the coexression correlation after shuffling
    corrShuffled <- c(corrShuffled,  stats::cor(gene1Exp, gene2Exp, use = "pairwise.complete.obs"))
  }
  return(corrShuffled)
}

#' \code{pairExpHist} Generate distribution of coexpression correlations between
#'                    randomly chosen gene pairs, compared to the distribution of
#'                    shuffled expression values for each pair of significantly
#'                    differetially coexpresssed genes
#'
#' \code{pairExpHist} Generate distribution of coexpression correlations between
#'                    randomly chosen gene pairs, compared to the distribution of
#'                    shuffled expression values for significantly differentially
#'                    coexpressed genes. In addition, this function displays the
#'                    distribution in a histogram for visualization and interpretation
#'
#' @param expProfile (data frame)  gene pairs with co expression correlation represented
#'                                 in a data frame
#'
#' @return (data frame) If there are significantly coexpressed gene pairs,
#'                      extract the expression profile data frame with two genes
#'                      that are significantly differentially coexpressed, their
#'                      correlation and p-value for their coexpression. Otherwise,
#'                      return the above information for all the gene pairs
#'
#' @importFrom stats cor
#' @import graphics
#'
#' @examples
#' \dontrun{
#' complex <- "MDC1-MRN-ATM-FANCD2 complex"
#' complexExpHeatmap(complex)
#' expProfile <- getCoExpProfile(complex)
#' pairExpHist(expProfile)
#' }
#' @export
#'
pairExpHist <- function(expProfile){
  N <- 100000
  myQNXPProcesed <- processmyQNXP()
  simulationCorr <- c()
  # Randomly generate 100000 gene pairs and their coexpression correlation
  for (i in 1:N){
    iGene1 <- sample(1:nrow(myQNXPProcesed), 1)
    iGene2 <- sample(1:nrow(myQNXPProcesed), 1)
    simulationCorr <- c(simulationCorr,  stats::cor(myQNXPProcesed[iGene1, ], myQNXPProcesed[iGene2, ], use = "pairwise.complete.obs"))
  }
  # calculate the p-Values of the correlation of the random gene pair distribution
  pValues <- apply(expProfile[3], 2, function(x) getpValue(simulationCorr, x))
  expProfile$pValues <- pValues
  colnames(expProfile)[4] <- "pValues"
  # record the gene pairs with a p-Values < 0.05
  iSigGenePairs <- which(pValues < 0.05)
  print((length(iSigGenePairs) != 0))
  if (length(iSigGenePairs) != 0){
    # generate histogram for significant gene pairs for visualization
    for (i in iSigGenePairs){
      print((length(iSigGenePairs) != 0))
      print(paste("genes", expProfile$Gene1[i], expProfile$Gene2[i]))
      genePairsShuffle <- genePairsShuffle(expProfile$Gene1[i], expProfile$Gene2[i])
      graphics::hist(genePairsShuffle, col="light pink", ylim = c(0, N/4), xlim = c(-1, 1),
           main = paste(expProfile$Gene1[i], " ", expProfile$Gene2[i], "Expression Profile Correlations"),
           ylab = "Density", xlab = "Coefficient of Correlation [1, -1]")
      graphics::hist(simulationCorr, col= "light blue", add=T)
      graphics::legend("topright", c("random", "genePair"), col=c("light blue", "light pink"), lwd=10)
    }
    # return information for significantly coexpressed gene pairs
    return(expProfile[iSigGenePairs, ])
  }else{
    print("No significant outlier in this complex")
    # if there is no significant gene pairs, return the who expression profile
    return(expProfile)
  }
}

# [END]

