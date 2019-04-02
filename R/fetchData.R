# fetchData.R

#' \code{fetchData} read a prepared dataset from a repository.
#'
#' \code{fetchData} loads and returns datasets from a data repository. When
#'                  invoked without parameters, a list of available datasets is
#'                  printed. When invoked with a parameter \code{set}, the
#'                  intended dataset is fetched with \code{\link[base]{readRDS}}
#'                  and returned.
#'
#' @section Details:
#' \itemize{
#'   \item For the \code{HGNCreference} dataset, see ...
#'   \item For the GTRD derived transcription factor datasets
#'           \code{GTRDgeneTFs} and \code{GTRDTFgenes}, see ...
#'   \item For the BioGRID protein-protein interaction dataset, see details in
#'             the \href{../doc/dataDetails-BioGRID.html}{BioGRID data vignette}
#'            (or load the vignette with \code{vignette("dataDetails-BioGRID",
#'            package = "BCB420.2019.ESA")}).
#'   \item For the  STRING database derived datasets \code{STRINGedges0.9},
#'            \code{STRINGedges0.8}, and \code{STRINGactions}, see ...
#'   \item For the InterPro database datasets \code{genesIPR}, and
#'            \code{IPRgenes} dataset, see ...
#'   \item For the expression profiles dataset \code{GEOprofiles}, see ...
#'   \item For the systems database dataset \code{SysDB}, see ...
#'   \item The Reactome dataset \code{ReactomeSym} is currently undocumented.
#' }
#'
#'
#' @param set (character)  name of the requested dataset
#'
#' @return a dataset, usually either a \code{list} or a \code{data frame}, or
#'         \code{NULL} (invisibly) if no dataset was specified.
#' @examples
#' fetchData()                         # prints available datasets
#' HGNC <- fetchData("HGNCreference")  # assigns the HGNC reference dataset
#' @export

fetchData <- function(set) {
  # ...
  baseURL <- "http://steipe.biochemistry.utoronto.ca/abc/assets/"

  availDat <- data.frame(sets = c("HGNCreference",
                                  "GTRDgeneTFs",
                                  "GTRDTFgenes",
                                  "BioGRID",
                                  "STRINGedges0.9",
                                  "STRINGedges0.8",
                                  "STRINGactions",
                                  "genesIPR",
                                  "IPRgenes",
                                  "GEOprofiles",
                                  "SysDB",
                                  "ReactomeSym"),
                         desc = c("HGNC symbols and crossreferences",
                                  "Genes and their TFs by ChIP-seq from GTRD",
                                  "TFs and their genes by ChIP-seq from GTRD",
                                  "BioGRID physical and genetic interaction data",
                                  "STRING edges with confidence p >= 0.9 mapped to HGNC symbols",
                                  "STRING edges with confidence p >= 0.8 mapped to HGNC symbols",
                                  "STRING actions mapped to HGNC symbols",
                                  "genes and the InterPro domains they contain",
                                  "InterPro domains and the genes they are found in",
                                  "GEO expression profiles, quantile normalized",
                                  "A systems database",
                                  "A tibble of Reactome IDs and HGNC symbols"),
                         FN   = c("HGNCreference.rds",
                                  "GTRDgeneTFs-2019-03-13.rds",
                                  "GTRDTFgenes-2019-03-13.rds",
                                  "BioGRID.3.5.170.rds",
                                  "STRINGedges-p0.9-2019-03-14.rds",
                                  "STRINGedges-p0.8-2019-03-26.rds",
                                  "STRINGactions-2019-03-26.rds",
                                  "genesIPR-V.73.rds",
                                  "IPRgenes-V.73.rds",
                                  "GEO-QN-profile-2019-03-24.rds",
                                  "systemsDB-2018-03-28.rds",
                                  "ReactomeSym.rds"),
                         stringsAsFactors = FALSE)
  rownames(availDat) <- availDat$sets

  if (missing(set)) {
    print(availDat[ , c("sets", "desc")])
  } else {
    x <- readRDS(url(paste0(baseURL, availDat[set, "FN"])))
    return(x)
  }
  return(invisible(NULL))
}

# [END]
