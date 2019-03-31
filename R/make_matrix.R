# make_matrix.R
#' Make Similarity Matrix
#'
#' \code{make_matrix.R} Generate a pairwise similarity matrix for a set of genes.
#'
#' The the similarity metric is determined by the similarity function that is passed as a parameter.
#' @section <title>: Additional explanation.
#'
#' @param dist_fn Function to be used to calculate individual pairwise similarities.
#' @param genes A vector of HGNC symbols for which to calculate the similarity of all possible pairwise combinaitons.
#' @param data_source An object that will be passed to the similarity function in order to calculate similarity of the 2 genes.
#' @return similarity matrix, containing numerics values between 0 and 1
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @seealso \code{\link{net_sim}} \code{\link{tf_sim}} \code{\link{expr_sim}} Possible similarity functions. See thier documentation to determine which data source needs to be provided for each.
#'
#' @examples
#' # Calculate a similarity matrix for BRCA1 and BRCA2 based on expression similarity
#' make_matrix(c("BRCA1", "BRCA2"), dist_fn = expr_sim, data_source = fetchData("GEOprofiles"))
#' # BRCA1      BRCA2
#' # BRCA1 1.00000000 0.05755273
#' # BRCA2 0.05755273 1.00000000
#'
#'
#' @export

make_matrix <- function(genes, dist_fn, data_source) {
  # initialize matrix
  l <- length(genes)
  dist_mat <- matrix(nrow = l, ncol = l)
  rownames(dist_mat) <- genes
  colnames(dist_mat) <- genes
  # populate the similarity matrix using the similarity function
  pb <- txtProgressBar(min = 0, max = length(genes), initial = 0, char = "=")
  cat("Building similarity matrix...\n")
  for (i in seq_along(genes)) {
    setTxtProgressBar(pb, i)
    for (j in seq_along(genes)) {
      dist_mat[i, j] <- dist_fn(genes[i], genes[j], data_source)
    }
  }
  return(dist_mat)
}

# [END]
