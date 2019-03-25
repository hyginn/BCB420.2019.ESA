#' Function semanticSimilarity().
# Function stub taken from https://github.com/hyginn/
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


#' \code{semanticSimilarity} outputs semantic similarity matrix of enriched
#' GO terms for a systems gene set
#'
#'
#' @param sys System Code
#' @return (matrix) a matrix of enriched GO terms with semantic similarity values
#'
#' @examples
#' semanticSimilarity("PHALY")
#'
#' @export
semanticSimilarity <- function(sys) {

  #Load the HGNC data
  myURL <- paste0("https://github.com/hyginn/",
                  "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))

  #org.Hs.eg.db uses entrez gene identifiers - need to annotate the geneset into entrez ID's with HGNC
  geneset <- fetchComponents(sys)
  geneset <- (HGNC$sym %in% geneset)
  genes <- HGNC[geneset,]$GeneID
  genes <- as.character(genes)

  #finding GO term enrichment using our geneset
  GOSim::setEvidenceLevel(evidences = "all",
                          organism = org.Hs.egORGANISM,
                          gomap = org.Hs.egGO)

  allGenes <- AnnotationDbi::keys(org.Hs.eg.db)
  GOSim::setOntology("BP", loadIC = FALSE)

  #enriched GO Terms from our sample gene set
  myEnr <- GOSim::GOenrichment(genes, allGenes)

  #semantic similarity of enriched GO terms as a matrix
  #(note that this will take really long accoring to the size of the data)
  semSim <- GOSemSim::mclusterSim(myEnr$genes,semData=hsGO, measure="Wang", combine="BMA")

  #returns of semantic similarity matrix of the enriched GO terms from our input geneset
  return(semSim)
}
