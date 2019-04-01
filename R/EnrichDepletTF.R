# EnrichDepletTF.R

# =====  1  Functions that loads the required datasets  ================================================
# ====== 1.1  Load expression-profiles  =============================

#' sysExProf.
#'
#' \code{sysExProf} Loads expression-profiles of chosen genes
#'
#' @section Details: The expression profiles were compiled by Boris Steipe from 52
#' microarray experiments downloaded from GEO, and quantile normalized. The function
#' that fetches the database was also written by Dr. Boris Steipe
#'\href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R}{fetchData Github}
#' The code that deletes experiments with na values was written by Boris Steipe
#'\href{ https://github.com/hyginn/BCB420.2019.ESA/blob/master/README.md}{removeNAs Github}
#'
#'
#' @param GeneSym (char vector)  A character vector of length > 0L of HGNC symbols.
#' @return  (matrix)             A matrix with rows dimentions of the GeneSym vector
#'                               length and 52 columns corresponding to the GEO experiments.
#'                               Retrieves the expression profiles of genes that have data
#'                               for all 52 GEO experiments. The rows names are the HGNC
#'                               sybmols and the column names are the experiments' names.
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' # Expression profiles for "BRCA1" and "AR" genes.
#' sysExProf(GeneSym = c("BRCA1", "AR"))
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
#' @section Details: The fetched dataset was collected and composed by Boris Steipe
#'\href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R{fetchData}.
#' The GTRDgeneTFs dataset lists which TFs binds the upstream regulatory regions
#' of each specified gene, as observed by ChIP-seq experiments to compose the GTRD
#' database. The loading function was written by Boris Steipe
#'\href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R{fetchData}.
#'
#' @param GeneSym (char vector)  A character vector of length > 0L of HGNC symbols.
#' @return (list)                A list of HGNC-symbol-named character vectors of
#'                               the transcription factors that bind each gene.
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' # Fetching the vector of the transcription factors that bind the gene BCL3.
#' sysGTRDgenes(GeneSym = c("BCL3"))
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
#' @section Details: The fetched dataset was collected and composed by Boris Steipe
#'\href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R{fetchData}.
#' The GTRDTFgenes dataset lists which genes bind each specified TF, as observed by
#' ChIP-seq experiments that compose the GTRD database. The loading function was written
#' by Boris Steipe
#'\href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R{fetchData}.
#'
#' @param TfSym (char vector)   A character vector of length > 0L of unique TFs
#'                              (HGNC symbols) that was extracted from the output
#'                              of sysGTRDgenes() for a certain set of genes.
#' @param GeneSym (char vector) A character vector of length > 0L of the unique genes
#'                              (HGNC symbols) that the function sysGTRDgenes()
#'                              retrieved for a certain set of genes (the same certain
#'                              set of genes that was used to produce TfSym).
#' @return (list)               A list of transcription-factor-named character vectors of
#'                              the genes that bind each transcription factor.
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' # Fetching the vector of the genes that bind the transcription factor AR.
#' GTRDgenes <- sysGTRDgenes(GeneSym = c("LAMP1", "HDAC6"))
#' GeneSym <- names(GTRDgenes) #The specific set of genes retrieved
#' TfSym <- unique(unlist(GTRDgenes, use.names = FALSE)) #TFs that bind the specific set of genes
#' sysGTRDtf(TfSym = TfSym, GeneSym = GeneSym)
#'

sysGTRDtf <- function(TfSym, GeneSym){

  GTRDtf <- fetchData("GTRDTFgenes")
  GTRDtf <- GTRDtf[which(names(GTRDtf) %in% TfSym)] #keeps only TFs that appear
                                                    #in TfSym (i.e., TF that binds
                                                    #the chosen genes.
  GTRDtf <- lapply(GTRDtf, function(x) x[which(x %in% GeneSym)]) #Keep only genes
                                                                 #that are in GeneSym

  #Check:
  if (identical(sort(unique(unlist(GTRDtf, use.names = FALSE))) ,sort(GeneSym))){
    return(GTRDtf)
  }else{
    return(print("Error: input genes are not equal to output genes"))
  }
}


# ===== 2 Finding positive correlated genes =============================

#' sysUpCorGenes.
#'
#' \code{sysUpCorGenes} Finding pairs of genes that their expression profiles are the
#'                      up most positively correlated.
#'
#' @section Details: A function that produces a character vector of genes (HGNC symbols)
#'                   of a set of unique genes that scored the up most correlation
#'                   values for the expression profiles.
#'
#' @param exProfS (matrix)     A matrix with 52 columns that corresponds to 52 microarray
#'                             experiments of expression profiles downloaded from GEO and
#'                             quantile normalized and fetched by Dr. Boris Steipe
#'                             \href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R}{fetchData Github}
#'                             The matrix is produced by the function sysExProf() and the
#'                             row dimentions are set by the length of the vector of the
#'                             queried genes that was fed into sysExProf().
#' @param nUpMost (numeric / boolean) Either "FALSE" to returns all significantly
#'                                    correlated genes, a numeric or vector of 1L length,
#'                                    describing the amount of upmost correlated pairs of
#'                                    genes to return by the function. The default is set
#'                                    to FALSE.
#' @param direction (char)     A character vector of 1L length of either "Positive" or
#'                             "Negative". "Positive" returns pairs of genes that are
#'                             signigicantly positively correlated. "Negative" returns
#'                             pairs of genes that are signigicantly negatively correlated.
#'                             The default is set to "Positive".
#' @param multipleTests (char) A character vector of 1L length of either "Bonferroni" or
#'                             "BH" (Benjamini-Hochberg), to set the multiple tests
#'                             correction approach. The default is set to "Bonferroni".
#' @param alpha (numeric)      a numeric or vector of 1L length, to set the value of alpha
#'                             for deciding on the cutoff for significant results. The
#'                             Default value is 0.05.
#' @return (char vector)       A character vector of unique correlated of genes (HGNC
#'                             symbols). Although the number of pairs of genes that
#'                             constitutes the vector is set by the paramenter "nUpMost",
#'                             the length of the vector may vary and won't always be
#'                             2 * "nUpMost", because duplications are deleted so that each
#'                             gene appears only once.
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{cor.test}} Test for association between paired samples, using one
#'                                 of Pearson's product moment correlation coefficient,
#'                                 Kendall's tau or Spearman's rho.
#'
#' @examples
#' # Vector of the 12 up most positively correlated pairs of genes of "SLIGR"
#' # system.
#' GeneSym <- SyDBgetSysSymbols(fetchData("SysDB"),"SLIGR")[[1]]
#' exProfs <- sysExProf(GeneSym = GeneSym)
#' sysUpCorGenes(exProfS = exProfs, nUpMost = 12)
#'

sysCorGenes <- function(exProfS,
                        nUpMost = FALSE,
                        direction = "Positive",
                        multipleTests = "Bonferroni",
                        alpha = 0.05){

  GeneSym <- row.names(exProfS)
  nPairs <- (length(GeneSym) * (length(GeneSym) - 1)) / 2 #number of pairs of genes
  expCor <- as.data.frame(matrix(nrow = nPairs, ncol = 6, data = c("")),
                          stringsAsFactors = FALSE)
  colnames(expCor) <- c("Gene1", "Gene2", "Correlation", "Pvalue", "BH", "Bonferroni")
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
  } else if (direction == "Negative") { #Choosing negatively correlated genes
    signifCor <- signifCor[signifCor$Correlation < 0, ]
  }

  if (nUpMost == FALSE){ #returning all the genes that are correlated
    MostCor <- signifCor[order(signifCor$Correlation, decreasing = TRUE), ]
  } else { #returning only the "nUpMost" correlated genes
  MostCor <-
    signifCor[order(signifCor$Correlation, decreasing = TRUE), ][1 : nUpMost, ]
  }

  return(unique(c(MostCor$Gene1, MostCor$Gene2)))
}

# ===== 3 Enrichement and depletion of TF binding sites =============================

#' EnrichDepletTF.
#'
#' \code{EnrichDepletTF} Calculating the p values for enrichment and depletion of
#'                       transcription factor binding sites within a system.
#'
#'
#' @section Details: The systems were curated as part of BCB420H1 winter 2019 course. The
#'                   TF binding datasets were curated by the GTRD database and edited
#'                   by Boris Steipe. The expression profiles database was obtained from
#'                   GEO and edited by Boris Steipe.
#'\href{https://github.com/hyginn/BCB420.2019.ESA/blob/master/R/fetchData.R}{fetchData Github}
#'                   Enrichment and depletion calculation were inspired by Boris Steipe's
#'                   guidance and by "Pathway Guide"
#'\href{https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/}{Fishers Exact Test}
#'
#' @param sys (char)           A character vector of 1L length of 5 uppercase lettered word
#'                             that corresponds to a biological system curated during the
#'                             BCB420H1 winter 2019 course. The list of available systems
#'                             can be retrieved by:
#'                             names(SyDBgetRootSysIDs(fetchData("SysDB")))
#' @param nUpMost (numeric / boolean) Either "FALSE" or numeric. Use "FALSE" for calculating
#'                                    one sided enrichment and one sided depletion for all
#'                                    significantly correlated genes. Use a a number (vector
#'                                    of 1L length) to choose an amount of upmost
#'                                    significantly correlated pairs of genes to calculate
#'                                    their enrichment / depletion. The default is "FALSE".
#' @param CorDirection (char)  Either "Positive" or "Negative". "Positive" uses pairs of
#'                             genes that are significantly positively correlated to
#'                             calculate their enrichment / depletion of TF binding.
#'                             "Negative" uses pairs of genes that are significantly
#'                             negatively correlated for the enrichment / depletion
#'                             analysis. The default is "Positive".
#' @param multipleTests (char) Either "Bonferroni" or "BH" (Benjamini-Hochberg). Sets the
#'                             multiple tests correction mathematical approach for the
#'                             calculation of correlation. The default is set to
#'                             "Bonferroni". Note: this parameter does not influence the
#'                             correction approach for the depletion / enrichment
#'                             calculations, for the returned argument presents both
#'                             approaches as default.
#' @param alpha (numeric)      A numeric vector of 1L length. Sets the value of alpha
#'                             as the cutoff for significant results. The same alpha is
#'                             used for both choosing significantly correlated genes and
#'                             for significantly enrichement / depletetion. The default
#'                             value is 0.05.
#' @return (data frame)        A data frame with 8 columns and rows number that is equal
#'                             to the amnout of unique TF that can bind the genes of the
#'                             chosen system. For each TF, provides: Enrichment magnitude
#'                             among the systems' correlated genes; P value for one sided
#'                             fishers exact test for enrichment; Depletion magnitude
#'                             among the systems' correlated genes; P value for one sided
#'                             fishers exact test for depletion; The BH cutoff according
#'                             to the p value for enrichment; The BH cutoff according to
#'                             the p value for depletion; The Bonferroni cutoff. Returns
#'                             all values for both significantly enriched / depleted TFs
#'                             and not-enriched / not-depleted TFs.
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0002-9478-5974}{Sapir Labes} (aut)
#'
#' @seealso \code{\link{fisher.test}} Performs Fisher's exact test for testing the null of
#'                                    independence of rows and columns in a contingency table with fixed marginals.
#'
#' @examples
#' # Creating the enrichment and depletion p values for TF binding by the positively
#' # correlated genes of HVGCR system.
#' EnrichDepletTF(sys = "HVGCR")
#'
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

  if (sum( ! is.na(upCorGenes)) >= 3) {
    #calculating enrichment and depletion:
    corEnrichDep <- as.data.frame(matrix(data = NA, nrow = length(TfSym), ncol = 8),
                                   stringsAsFactors = F)
    colnames(corEnrichDep) <- c("TF", "Enrichment", "Enrichment_P_value",
                                "Depletion", "Depletion_P_value",
                                "BH_enrichment","BH_depletion",
                                "Bonferroni")
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

  } else {
    return (print("Error: less then 3 significantly correlated genes"))
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
#if (identical(sort(unique(unlist(GTRDtf, use.names = FALSE))) ,sort(GeneSym))){
#  return(GTRDtf)
#}else{
#  return(print("Error: input genes are not equal to output genes"))
#}



# [END]
