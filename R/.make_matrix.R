# make_matrix.R

# Helper function for clusterSystems
# probably not useful elsewhere

#' Make Distance Matrix
#'
#' \code{make_matrix.R} Generate a pairwise distance matrix for a set of genes.
#'
#' The the distance metric is determined by the distance function that is passed as a parameter.
#' The main use of this function is as a helper function for \code{\link{clusterSystems}}.
#'
#' @param dist_fn Function to be used to calculate individual pairwise distances.
#' @param genes A vector of HGNC symbols for which to calculate the similarity of all possible pairwise combinaitons.
#' @param data_source An object that will be passed to the similarity function in order to calculate similarity of the 2 genes.
#' @return Object of class distance
#'

#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @seealso \code{\link{jaccard_dist}}, \code{\link{tf_dist}} and \code{\link{expr_dist}} are possible similarity functions that can be used as input. See thier documentation to determine which data source needs to be provided for each.
#'
#' @examples
#' # Calculate a similarity matrix for BRCA1 and BRCA2 and PTEN based on expression similarity
#' data_source <- fetchData("GEOprofiles")
#' matrix <- make_matrix(c('BRCA1', 'BRCA2', 'PTEN'), expr_dist, data_source)
#' #           BRCA1     BRCA2
#' # BRCA2 0.9424473
#' # PTEN  0.7862693 0.3481342
#'
#'


make_matrix <- function(genes, dist_fn, data_source) {
  # initialize matrix
  l <- length(genes)
  dist_mat <- matrix(nrow = l, ncol = l)
  rownames(dist_mat) <- genes
  colnames(dist_mat) <- genes
  # populate the similarity matrix using the similarity function
  pb <- utils::txtProgressBar(min = 0, max = length(genes), initial = 0, char = "=")
  cat("\n\nBuilding distance matrix...\n")
  for (i in seq_along(genes)) {
    utils::setTxtProgressBar(pb, i)
    for (j in seq_along(genes)) {
      dist_mat[i, j] <- dist_fn(genes[i], genes[j], data_source)
    }
  }
  dist_mat <- stats::as.dist(dist_mat) # distance class object
  return(dist_mat)
}

# [END]
