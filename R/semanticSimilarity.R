# semanticSimilarity.R
#'
#' \code{semanticSimilarity} gets semantic similarity matrix of enriched
#' GO terms for a systems gene set
#'
#'
#' @param sys System Code
#' @return (matrix) a matrix of enriched GO terms with semantic similarity values
#'
#' @author \href{https://orcid.org/0000-0002-1134-6758}{Boris Steipe} (aut)
#'
#' @import BiocManager
#' @import BiocCheck
#' @import AnnotationDbi
#' @import GOSim
#' @import GOSemSim
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

# [END]
