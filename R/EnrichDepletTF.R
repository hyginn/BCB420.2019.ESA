# EnrichDepletTF.R

# =====  1  Functions that loads the required datasets  ================================================
# ====== 1.1  Load expression-profiles  =============================

#' sysExProf.
#'
#' \code{sysExProf} Loads expression-profiles of chosen genes
#'
#' @section Details: This is a helper function. The expression profiles were compiled
#' from microarray experiments downloaded from GEO, and quantile normalized. For the
#' function that fetches the database, see \code{\link{fetchData}}.
#'
#' @param GeneSym (character)   A vector of length > 0L of HGNC symbols.
#' @return  (numeric)           A matrix with row dimensions that are equal to the
#'                              number of genes of the GeneSym vector that has
#'                              expression data for all of the experiments.
#'                              The column dimensions corresponds to the number of
#'                              microarray experiments in the expression profile
#'                              dataset. The matrix contains the expression profiles
#'                              of genes that have data for all experiments in the
#'                              expression profile source dataset. The rows names
#'                              are the HGNC sybmols and the column names are the
#'                              experiments' names.
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{fetchData}} Fetches the requested dataset.
#'
#' @examples
#' # Expression profiles for "BRCA1" and "AR" genes.
#' # sysExProf(GeneSym = c("BRCA1", "AR"))
#'

sysExProf <- function(GeneSym) {

  # Load the expression profiles:
  myQNXP <- fetchData("GEOprofiles")

  #Choosing genes that have data for all experiments
  GeneSym <- GeneSym[( ! is.na(rowMeans(myQNXP[GeneSym,], na.rm = FALSE)))]
  myQNXP <- myQNXP[GeneSym, ]

  return(myQNXP)
}

# ====== 1.2  Loading the GTRDgeneTFs data for chosen genes =============================

#' sysGTRDgenes.
#'
#' \code{sysGTRDgenes} Fetch GTRD data for the specified genes.
#'
#' @section Details: This is a helper function. The fetched dataset was composed from the
#' GTRDgeneTFs dataset which lists which TFs bind the upstream regulatory regions of each
#' specified gene, as observed by ChIP-seq experiments. For more details about the
#' function that fetches the GRTDgeneTFs database, see \code{\link{fetchData}}.
#'
#' @param GeneSym (character)  A vector of length > 0L of HGNC symbols.
#' @return (list)              A list of HGNC-symbol-named character vectors of
#'                             the transcription factors that bind each gene.
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{fetchData}} Fetches the requested dataset.
#'
#' @examples
#' # Fetching the vector of the transcription factors that bind the gene BCL3.
#' # sysGTRDgenes(GeneSym = c("BCL3"))
#'

sysGTRDgenes <- function(GeneSym){

  GTRDgenes <- fetchData("GTRDgeneTFs")
  GTRDgenes <- GTRDgenes[which(names(GTRDgenes) %in% GeneSym)] #Keeping the data for
                                                               #the specified genes
  return(GTRDgenes)
}

# ====== 1.3  Loading the GTRDTFgenes data for chosen genes and TF =============================

#' sysGTRDtf.
#'
#' \code{sysGTRDtf} Fetch GTRD data for the specified Transcription Factors and genes.
#'
#' @section Details: This is a helper function. The fetched dataset was collected and
#' composed from the GTRD database. The GTRDTFgenes dataset lists which  genes bind
#' each specified TF, as observed by ChIP-seq experiments. For more details about the
#' function that fetches the GTRDTFgenes database, see \code{\link{fetchData}}. This
#' function requires an input that can be produced by the helper function sysGTRDgenes.
#'
#' @param TfSym (character)     A vector of length > 0L of unique TFs (HGNC symbols)
#'                              that was extracted from the output of sysGTRDgenes()
#'                              for a certain set of genes.
#' @param GeneSym (character)   A vector of length > 0L of the unique genes (HGNC symbols)
#'                              that the function sysGTRDgenes() retrieved for an input
#'                              of a certain set of genes (the same certain set of genes
#'                              that was used to produce TfSym).
#' @return (list)               A list of transcription-factor-named character vectors of
#'                              the genes that bind each transcription factor.
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{fetchData}} Fetches the requested dataset.
#'
#' @examples
#' \dontrun{
#' # Fetching the vector of the genes that bind the transcription factor AR.
#' GTRDgenes <- sysGTRDgenes(GeneSym = c("LAMP1", "HDAC6"))
#' mySymbols <- names(GTRDgenes) #The specific set of genes retrieved
#' myTfSym <- unique(unlist(GTRDgenes, use.names = FALSE)) #TFs that bind the specific set of genes
#' sysGTRDtf(TfSym = myTfSym, GeneSym = mySymbols)
#' }

sysGTRDtf <- function(TfSym, GeneSym){

  GTRDtf <- fetchData("GTRDTFgenes")
  GTRDtf <- GTRDtf[which(names(GTRDtf) %in% TfSym)] #keeps only TFs that appear
                                                    #in TfSym (i.e., TF that binds
                                                    #the chosen genes.
  GTRDtf <- lapply(GTRDtf, function(x) x[which(x %in% GeneSym)]) #Keep only genes
                                                                 #that are in GeneSym

  #Check:
  if (! identical(sort(unique(unlist(GTRDtf, use.names = FALSE))) ,sort(GeneSym))){
    stop(sprintf("%s is not identical to the genes that bind 'TfSym'.", "'GeneSym'"))
  }
    return(GTRDtf)
}


# ===== 2 Finding positively or negatively correlated genes =============================

#' sysUpCorGenes.
#'
#' \code{sysUpCorGenes} Finding pairs of genes that their expression profiles are the
#'                      most positively or negatively correlated.
#'
#' @section Details: A helper function that produces a character vector of unique genes
#'                   (HGNC symbols) that belongs to pairs of genes that their expression
#'                   profiles were shown to be either positively or negatively correlated.
#'                   It is also possible to specify the number the of most highly
#'                   negatively or positively correlated pairs of genes to retrieve.
#'
#' @param exProfS (numeric)    A matrix of expression profiles, with column dimensions
#'                             equal to the number of microarray experiments, and row
#'                             dimensions that are equal to the number of queried genes.
#'                             All genes must have information for all experiments (i.e.,
#'                             no NAs are allowed). The matrix can be produced  by the
#'                             helper function sysExProf().
#' @param nUpMost (numeric|logical) Either FALSE or a numeric vector of 1L length. FALSE
#'                                  returns all significantly correlated genes. A numeric
#'                                  vector sets the number of most highly correlated
#'                                  pairs of genes to return by the function. The default
#'                                  is FALSE.
#' @param direction (character) A vector of 1L length of either "Positive" or "Negative".
#'                              "Positive" returns pairs of genes that are significantly
#'                              positively correlated. "Negative" returns pairs of genes
#'                              that are significantly negatively correlated. The
#'                              default is set to "Positive".
#' @param multipleTests (character) A vector of 1L length of either "Bonferroni" or
#'                                  "BH" (Benjamini-Hochberg). Sets the multiple tests
#'                                  correction approach. The default is "Bonferroni".
#' @param alpha (numeric)      A vector of 1L length. Sets the value of alpha for
#'                             deciding the cutoff for significant results. The Default
#'                             is 0.05.
#' @return (character)         A vector of unique genes (HGNC symbols) that were found to
#'                             be significantly correlated, according to the conditions
#'                             set by the parameters. For a numeric value of nUpMost,
#'                             although the number of pairs of genes was set, the number
#'                             of the returned genes is not always equal to 2 * "nUpMost",
#'                             because some of the pairs may contain the same genes,
#'                             while the function returns only a vector of unique genes.
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{cor.test}} Test for association between paired samples, using one
#'                                 of Pearson's product moment correlation coefficient,
#'                                 Kendall's tau or Spearman's rho.
#'
#' @examples
#' \dontrun{
#' # Vector of the 12 up most positively correlated pairs of genes of "SLIGR"
#' # system.
#' # mySymbols <- SyDBgetSysSymbols(fetchData("SysDB"),"SLIGR")[[1]]
#' # exProfS <- matrix(data = sample(1:200, 4 * length(mySymbols), replace = TRUE),
#' #                  ncol = 4) #generating synthetic expression data
#' # row.names(exProfS) <- mySymbols
#' # sysUpCorGenes(exProfS = exProfS, nUpMost = 12)
#' }

sysCorGenes <- function(exProfS,
                        nUpMost = FALSE,
                        direction = "Positive",
                        multipleTests = "Bonferroni",
                        alpha = 0.05){

  #checks
  if (sum(is.na(exProfS)) != 0){ #No NAs in the exProfS matrix
    stop(sprintf("%s contains NAs.", "'exProfS'"))
  }
  if (!((direction == "Positive") || (direction == "Negative"))){
    stop(sprintf("%s must be either 'Positive' or 'Negative'.", "'direction'"))
  }
  if (!((multipleTests == "Bonferroni") || (multipleTests == "BH"))){
    stop(sprintf("%s must be either 'Bonferroni' or 'BH'.", "'multipleTests'"))
  }

  GeneSym <- row.names(exProfS)
  nPairs <- (length(GeneSym) * (length(GeneSym) - 1)) / 2 #number of pairs of genes
  expCor <- data.frame(Gene1 = rep("", nPairs), Gene2 = rep("", nPairs),
                       Correlation = rep("", nPairs), Pvalue = rep("", nPairs),
                       BH = rep("", nPairs), Bonferroni = rep("", nPairs),
                       stringsAsFactors = FALSE)
  expCor$Gene1 <- unlist(lapply(GeneSym, #Adding the column of the first genes of the pairs.
                                function (x) rep(x = x,
                                                 (length(GeneSym) - which(GeneSym == x)))))
  expCor$Gene2 <- unlist(lapply(GeneSym[2 : length(GeneSym)], #Adding the seconde genes
                                function (x) GeneSym[(which(GeneSym == x)) : length(GeneSym)]))

  expCor$Correlation <- #Calculating the correlation between expression profiles.
    apply(expCor, 1, function(x) stats::cor.test(exProfS[x["Gene1"], ],
                                                 exProfS[x["Gene2"], ])$estimate)
  expCor$Pvalue <- #Calculating the P values.
    apply(expCor, 1, function(x) stats::cor.test(exProfS[x["Gene1"], ],
                                                 exProfS[x["Gene2"], ])$p.value)

  expCor$BH <- (order(expCor$Pvalue, decreasing = FALSE) / nPairs) * alpha #Cutoff value
                                                                      #for BH correction
  expCor$Bonferroni <- alpha / nPairs #Cutoff value for Bonferroni correction

  signifCor <-
    expCor[which(expCor$Pvalue < expCor[, multipleTests]), ] #Genes that are significantly
                                                            #correlated according to the
                                                           #selected multipleTests
                                                           #correction approach.
  if (direction == "Positive") { #Choosing positively correlated genes
    signifCor <- signifCor[signifCor$Correlation > 0, ]
    if (nUpMost == FALSE){ #returning all the genes that are correlated
      MostCor <- signifCor
    } else { #returning only the "nUpMost" correlated genes
      MostCor <-
        signifCor[order(signifCor$Correlation, decreasing = TRUE), ][1 : nUpMost, ]
    }
  } else if (direction == "Negative") { #Choosing negatively correlated genes
    signifCor <- signifCor[signifCor$Correlation < 0, ]
    if (nUpMost == FALSE){ #returning all the genes that are correlated
      MostCor <- signifCor
    } else { #returning only the "nUpMost" correlated genes
      MostCor <-
        signifCor[order(signifCor$Correlation, decreasing = FALSE), ][1 : nUpMost, ]
    }
  }
  return(unique(c(MostCor$Gene1, MostCor$Gene2)))
}

# ===== 3 Enrichment and depletion of TF binding sites =============================

#' EnrichDepletTF.
#'
#' \code{EnrichDepletTF} Calculating the p values for enrichment and depletion of
#'                       transcription factor binding sites within a system.
#'
#' @section Details: The systems were curated as part of BCB420H1 winter 2019 course. The
#'                   TF binding datasets were curated by the GTRD database. The
#'                   expression profiles database was obtained from GEO. For more
#'                   information about the datasets, see \code{\link{fetchData}}.
#'                   Enrichment and depletion calculation were inspired by "Pathway
#'                   Guide" \href{https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/}{Fishers Exact Test}
#'                   and by Devon Ryan \href{https://www.biostars.org/p/102946/}{Depletion}
#'
#' @param sys (character) A vector of 1L length, of 5 uppercase lettered word that
#'                        corresponds to a biological system curated during the
#'                        BCB420H1 winter 2019 course. The list of available systems
#'                        can be retrieved by:
#'                        names(SyDBgetRootSysIDs(fetchData("SysDB"))). For more
#'                        information, see \code{\link{fetchData}} and
#'                        \code{\link{SyDBgetRootSysIDs}}
#' @param nUpMost (numeric|logical) Either FALSE or a numeric vector of 1L length. FALSE
#'                                  returns the one sided enrichment / depletion scores
#'                                  and p values for all significantly correlated genes.
#'                                  A numeric vector sets the number of most highly
#'                                  correlated pairs of genes to return their one sided
#'                                  enrichment / depletion scores and p values.
#'                                  The default is FALSE.
#' @param CorDirection (character) A vector of 1L length of either "Positive" or
#'                                 "Negative". "Positive" calculates the enrichment /
#'                                 depletion of TF-binding among pairs of genes that are
#'                                 significantly positively correlated. "Negative"
#'                                 calculates the enrichment / depletion of TF-binding
#'                                 among pairs of genes that are significantly negatively
#'                                 correlated. The default is "Positive".
#' @param multipleTests (character) A vector of 1L length of either "Bonferroni" or
#'                                  "BH" (Benjamini-Hochberg). Sets the multiple tests
#'                                  correction approach for the calculation of
#'                                  correlation. The default is "Bonferroni". Note that
#'                                  this parameter does not influence the correction
#'                                  approach for the depletion / enrichment calculations,
#'                                  as the returned value contains both approaches as
#'                                  a default.
#' @param alpha (numeric)      A vector of 1L length. Sets the value of alpha for
#'                             deciding the cutoff for significant results. The same
#'                             value of alpha is used for both calculating the cutoff
#'                             p value for significantly correlated pairs of genes and
#'                             for calculating the cutoff p value for significantly
#'                             enriched / depleted TFs among the significantly correlated
#'                             pairs of genes. The Default is 0.05.
#' @return (data frame)        A data frame with 8 columns, and rows dimensions that are
#'                             equal to the number of unique TF that can bind the genes of
#'                             the chosen system. The returned data frame provides the
#'                             following information for each TF: Enrichment magnitude
#'                             (Odds Ratio) among the systems' correlated genes; P value
#'                             for one sided fisher's exact test for enrichment; Depletion
#'                             magnitude (Odds Ratio) among the systems' correlated genes;
#'                             P value for one sided fisher's exact test for depletion;
#'                             The BH corrected cutoff calculated for the p values of
#'                             the enrichment; The BH corrected cutoff calculated for
#'                             the p values of depletion; The Bonferroni corrected cutoff.
#'                             This data frame contains all values for the significantly
#'                             enriched / depleted TFs as well as for the not-enriched /
#'                             not-depleted TFs. If less then 3 correlated genes were fed
#'                             into the function, the function will return NaN and a
#'                             warning.
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{fisher.test}} Performs Fisher's exact test for testing the null of
#'                                    independence of rows and columns in a contingency table
#'                                    with fixed marginals.
#'
#' @examples
#' \dontrun{
#' # The enrichment and depletion values and p values of TF for the 10 most positively
#' # correlated genes of NLRIN system.
#' EnrichDepletNLRIN <- EnrichDepletTF(sys = "NLRIN", nUpMost = 10)
#' EnrichDepletNLRIN[EnrichDepletNLRIN$Enrichment_P_value < EnrichDepletNLRIN$BH_enrichment, ]
#' }
#' @export


EnrichDepletTF <- function(sys,
                       nUpMost = FALSE,
                       CorDirection = "Positive",
                       multipleTests = "Bonferroni",
                       alpha = 0.05){

  GeneSym <- SyDBgetSysSymbols(fetchData("SysDB"), sys)[[1]] #Fetching the genes
                                                             #of the system
  exProfS <- sysExProf(GeneSym = GeneSym) #Fetching the expression profiles
  GeneSym <- row.names(exProfS) #genes that have expression data
  GTRDgenes <- sysGTRDgenes(GeneSym = GeneSym) #Fetch GTRDgeneTFs data
  GeneSym <- names(GTRDgenes) #genes that have data on GTRDgeneTFs (= bind any TF)
  TfSym <- unique(unlist(GTRDgenes, use.names = FALSE)) #TFs that bind the the chosen
  #genes.
  GTRDtf <- sysGTRDtf(GeneSym = GeneSym, TfSym = TfSym) #Fetch GTRDTFgenes data for
                                                        #The system's TFs and genes.
  exProfS <- exProfS[GeneSym, ] #Keep expression profiles of genes that have GTRD data.
  upCorGenes <- sysCorGenes(exProfS = exProfS,    #the genes that were
                            nUpMost = nUpMost,    #Upmost correlated.
                            direction = CorDirection,
                            multipleTests = multipleTests,
                            alpha = alpha)

  if (sum( ! is.na(upCorGenes)) >= 3) { #3 or more correlated genes.
    #calculating enrichment and depletion:
    numRow <- length(TfSym)
    corEnrichDep<-  data.frame(TF = rep(NA, numRow),
                               Enrichment = rep(NA, numRow),
                               Enrichment_P_value = rep(NA, numRow),
                               Depletion = rep(NA, numRow),
                               Depletion_P_value = rep(NA, numRow),
                               BH_enrichment = rep(NA, numRow),
                               BH_depletion = rep(NA, numRow),
                               Bonferroni = rep(NA, numRow),
                               stringsAsFactors = F)
    corEnrichDep$TF <- TfSym

    for (i in seq_along(TfSym)){
      #Calculating enrichment and p value
      GenesBindTf <- GTRDtf[[TfSym[i]]]
      CorGenesBindTf <- upCorGenes[which(upCorGenes %in% GenesBindTf)]
      a <- length(CorGenesBindTf) #Correlates and binds TF
      b <- length(GenesBindTf) - a #not correlates and binds TF
      c <- length(upCorGenes) - a #Correlates and don't bind
      d <- (length(GeneSym) - length(upCorGenes)) - b #Not correlates and not binds

      tmpFisher <- stats::fisher.test(matrix(c(a, b, c, d), nrow = 2),
                                      alternative = "greater")
      corEnrichDep$Enrichment[i] <- tmpFisher$estimate
      corEnrichDep$Enrichment_P_value[i] <- tmpFisher$p.value

      #Calculating depletion. and p value
      tmpFisher <- stats::fisher.test(matrix(c(c, d, a, b), nrow = 2),
                                      alternative = "greater")
      corEnrichDep$Depletion[i] <- tmpFisher$estimate
      corEnrichDep$Depletion_P_value[i] <- tmpFisher$p.value
    }

    #Calculating the cutoffs for each multi-test correction approach
    nExperiments <- length(TfSym)
    corEnrichDep <- corEnrichDep[order(corEnrichDep$Enrichment_P_value,
                                         decreasing = FALSE),] #order according to
                                                               #increasing p values
    corEnrichDep$BH_enrichment <- (order(corEnrichDep$Enrichment_P_value,
                                          decreasing = FALSE) / nExperiments) * alpha
    corEnrichDep <- corEnrichDep[order(corEnrichDep$Depletion_P_value,
                                         decreasing = FALSE),] #order according to
                                                               #increasing P values
    corEnrichDep$BH_depletion <- (order(corEnrichDep$Depletion_P_value,
                                         decreasing = FALSE) / nExperiments) * alpha
    corEnrichDep$Bonferroni <- alpha / nExperiments #Bonferroni cutoff value.

    return(corEnrichDep)

  } else { #Less then 3 correlated genes.
    return(message(NaN))
  }
}


##Tests that I used during the process of writing the code:
##sysCorGenes()
##After creating the columns of the pairs of genes, checking that no pair is duplicated:
#if (sum(duplicated(expCor[,1:2])) != 0){
#  print("Error! Duplicated pairs of genes")
#}

##Testing whether the genes that appears in the output list of sysGTRDtf are
##the exact same genes that were given as the input in GeneSym
#if (! identical(sort(unique(unlist(GTRDtf, use.names = FALSE))) ,sort(GeneSym))){
#  stop(sprintf("%s is not identical to the genes that bind 'TfSym'.", "'GeneSym'"))
#}
#return(GTRDtf)
#}



# [END]
