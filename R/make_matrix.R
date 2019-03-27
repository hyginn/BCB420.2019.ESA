# make_matrix.R
#' Make Similarity Matrix
#'
#' \code{make_matrix.R} Generate a pairwise similarity matrix for a set of genes.
#'
#' The the similarity metric is determined by the similarity function that is passed as a parameter.
#' @section <title>: Additional explanation.
#'
#' @param sim_fn Function to be used to calculate individual pairwise similarities.
#' @param genes A vector of HGNC symbols for which to calculate the similarity of all possible pairwise combinaitons.
#' @return similarity matrix, containing numerics values between 0 and 1
#'
#' @family <optional description of family>
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @seealso \code{\link{<function>}} <describe related function>, ... .
#'
#' @examples
#' # <explain what the example does>
#'
#'
#' @export

make_matrix <- function(genes, sim_fn) {
  # initialize matrix
  l <- length(genes)
  sim_mat <- matrix(nrow = l, ncol = l)
  rownames(sim_mat) <- genes
  colnames(sim_mat) <- genes
  # populate the similarity matrix using the similarity function
  for (i in seq_along(genes)) {
    for (j in seq_along(genes)) {
      sim_mat[i, j] <- sim_fn(genes[i], genes[j])
    }
  }
  return(sim_mat)
}

# [END]
