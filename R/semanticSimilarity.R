# semanticSimilarity.R
#'
#' \code{semanticSimilarity} gets semantic similarity matrix of enriched
#' GO terms for a systems gene set
#'
#'
#' @param sys - System Code (e.g., PHALY)
#' @return (matrix) a matrix of enriched GO terms with semantic similarity values or null if system does not contain any genes (GOterm x GOterm format for the matrix)
#'
#' @author \href{https://orcid.org/0000-0003-4609-4965}{Cathy Cha} (aut)
#' @seealso Documentation of semantic similarity calculations and GO enrichment analysis \code{\link{GOSemSim}}
#'
#'
#' @examples
#' \dontrun{
#' semanticSimilarity("PHALY")
#' }
#' @export
semanticSimilarity <- function(sys) {

  #Load the HGNC data
  HGNC <- fetchData("HGNCreference")
  #imports the Genome wide annotation for Human
  #in separate helper function so that it is only imported when this function is used
  HuGenAnnotImport()

  #org.Hs.eg.db uses entrez gene identifiers - need to annotate the geneset into entrez ID's with HGNC
  geneset <- SyDBgetSysSymbols(HGNC, sys)
  geneset <- (HGNC$sym %in% geneset)
  genes <- HGNC[geneset,]$GeneID
  genes <- as.character(genes)

  if (length(genes) == 0) {
    return(NULL)
  }

  #finding GO term enrichment using our geneset
  GOSim::setEvidenceLevel(evidences = "all",
                          organism = "Homo sapiens",
                          gomap = org.Hs.egGO)

  allGenes <- AnnotationDbi::keys(org.Hs.eg.db)

  #enriched GO Terms from our sample gene set
  myEnr <- GOSim::GOenrichment(genes, allGenes)

  #semantic similarity of enriched GO terms as a matrix
  #(note that this will take really long accoring to the size of the data)
  semSim <- GOSemSim::mclusterSim(myEnr$genes,semData=hsGO, measure="Wang", combine="BMA")

  #returns of semantic similarity matrix of the enriched GO terms from our input geneset
  return(semSim)
}

# [END]
