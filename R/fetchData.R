# fetchData.R

#' \code{fetchData} read a prepared dataset from a repository.
#'
#' \code{fetchData} documentation forthcoming.
#'
#' @param set (character)  name of the requested dataset
#'
#' @return a dataset
#' @examples
#' fetchData()                         # prints avaialable datasets
#' HGNC <- fetchData("HGNCreference")  # assigns the HGNC reference dataset
#' @export

fetchData <- function(set) {
  # ...
  baseURL <- "http://steipe.biochemistry.utoronto.ca/abc/assets/"

  availDat <- data.frame(sets = c("HGNCreference",
                                  "GTRDgeneTFs",
                                  "GTRDTFgenes",
                                  "STRINGedges0.9",
                                  "STRINGedges0.8",
                                  "STRINGactions",
                                  "genesIPR",
                                  "IPRgenes",
                                  "GEOprofiles",
                                  "SysDB"),
                         desc = c("HGNC symbols and crossreferences",
                                  "Genes and their TFs by ChIP-seq from GTRD",
                                  "TFs and their genes by ChIP-seq from GTRD",
                                  "STRING edges with confidence p >= 0.9 mapped to HGNC symbols",
                                  "STRING edges with confidence p >= 0.8 mapped to HGNC symbols",
                                  "STRING actions mapped to HGNC symbols",
                                  "genes and the InterPro domains they contain",
                                  "InterPro domains and the genes they are found in",
                                  "GEO expression profiles, quantile normalized",
                                  "A systems database"),
                         FN   = c("HGNCreference.rds",
                                  "GTRDgeneTFs-2019-03-13.rds",
                                  "GTRDTFgenes-2019-03-13.rds",
                                  "STRINGedges-p0.9-2019-03-14.rds",
                                  "STRINGedges-p0.8-2019-03-26.rds",
                                  "STRINGactions-2019-03-26.rds",
                                  "genesIPR-V.73.rds",
                                  "IPRgenes-V.73.rds",
                                  "GEO-QN-profile-2019-03-24.rds",
                                  "systemsDB-2018-03-27.rds"),
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
