# domainEnrichment.R

#' Domain enrichment in a system.
#'
#' \code{domainEnrichment} Find the enriched domains in selected system and interpret enrichment by plot.
#'
#' Details.
#' @section <title>: Additional explanation.
#'
#' @param system The shortehand notation of selected system, five-letter code as discussed in class.
#' @return <description>.
#'
#' @author \href{https://orcid.org/0000-0001-8720-5874}{Fan Shen} (aut)
#'
#' @examples
#' system <- "PHALY"
#' domainEnrichment(system)
#'
#' @export
domainEnrichment <- function(system = "PHALY") {
  # 1. Fetch data from Dr. Steipe's BCB420.2019.ESA package
  # HGNC data
  myURL <- paste0("https://github.com/hyginn/",
    "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))

  # Genes and their domains
  genesIPR <- fetchData("genesIPR")

  # domains and their genes
  IPRgenes <- fetchData("IPRgenes")

  # Get the genes for the system
  genes <- fetchComponents(system)


  # 2. Find the enrichment
  # 2.1 store all the domains that in the system and their enrichment
  domains <- c()
  for (gene in genes) {
    domains <- c(genesIPR[gene][[1]], domains)
  }

  # 2.1 find the enrichemnt for each domain in the system

  # 2.2 fidn the enrichment for each domain in all genes


  # 3. Find the description for enriched domains
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  filters <- biomaRt::listFilters(myMart)
  filters[grep("HGNC", filters$description), ]
  attributes <- biomaRt::listAttributes(myMart)
  attributes[grep("Interpro", attributes$description), ]

  tmp <- biomaRt::getBM(filters = "hgnc_id",
    attributes = c("interpro_short_description","interpro_description"),
    values = genes,
    mart = myMart)

  # Multiple testing

  # plot

  # TODO: Wrap up
  # describe the example a little bit!
  # use devtools::document() to get the namespace and md file
  # Add pkg for description
}

# [END]
