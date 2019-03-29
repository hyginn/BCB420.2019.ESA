# From Dr. Steipe's BCB420.2019.ESA package
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
  }   else {
    s <- ""
  }
  return(s)
}


#' \code{seqSimVis()} Calculates sequence similarity and plots dot plots based on the
#' three highest similarity values to visualize
#'
#'
#' @param gene_of_int (char) gene of interest
#' @param system (char) system of interest

#' @examples
#' seqSimVis("BRCA1", "PHALY")
#'
#' @export
seqSimVis <- function(gene_of_int, system) {

  # set mart
  ensembl <- biomaRt::useMart("ensembl",
                        dataset = "hsapiens_gene_ensembl",
                        host    = "grch37.ensembl.org",
                        path    = "/biomart/martservice")

  # get sequence of gene of interest
  print(paste0("Getting sequence for gene of interest: ", gene_of_int))
  seq1 <- biomaRt::getSequence(id = gene_of_int,
                      type = "hgnc_symbol",
                      seqType = "peptide",
                      mart = ensembl)$peptide
  for (seq in seq1) {
    if (seq != "Sequence unavailable") {
      seq1 <- seq
      break
    }
  }

  # get genes of system
  sys_genes <- fetchComponents(system)

  vals <- c()
  sys_sequences <- c()
  # get similarity values of genes in system
  for (gene in sys_genes) {
    print(paste0("Getting sequence for gene: ", gene))
    # get sequence of gene
    seq2 <- biomaRt::getSequence(id = gene,
                        type = "hgnc_symbol",
                        seqType = "peptide",
                        mart = ensembl)$peptide

    for (seq in seq2) {
      if (seq != "Sequence unavailable") {
        seq2 <- seq
        break
      }
    }

    if (length(seq2) == 0) {
      sys_sequences <- c(sys_sequences, "")
    } else {
      sys_sequences <- c(sys_sequences, seq2)
    }

    # get similarity values
    plist <- list(seq1, seq2)
    if (length(seq2) != 0) {
      psimmat = protr::parSeqSim(plist, cores = 4, type = "local", submat = "BLOSUM62")
      print(paste0("Similarity value comparing against ", gene, ": ", psimmat[1,2]))
      vals <- c(vals, psimmat[1,2])
    } else {
      vals <- c(vals, 0)
      print(paste0("Similarity value comparing against ", gene, ": ", 0))
    }

    # because parSeqSim runs in parallel
    closeAllConnections()
  }

  par(mfrow=c(2,2))
  # get genes of top 3 values to plot
  gene_df <- data.frame(sys_genes, vals, sys_sequences, stringsAsFactors = F)
  gene_df <- gene_df[order(-vals),]

  print("Plotting...")
  # plot
  plot(vals, col=ifelse(vals %in% gene_df$vals[1:3], 'red', 'gray'), main="Similarity Values")

  # gene of interest char sequence
  seq1_chars <- unlist(strsplit(seq1,""))

  for (i in 1:3) {
    seq2_chars <- unlist(strsplit(gene_df$sys_sequences[i], ""))
    seqinr::dotPlot(seq1_chars, seq2_chars, xlab=gene_of_int, ylab=gene_df$sys_genes[i], main=paste0("Sim Value: ", gene_df$vals[i]))
  }

  print("..Done")

}


