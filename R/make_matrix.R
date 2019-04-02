# make_matrix.R

# Purpose: Generate a pairwise distance matrix for a set of genes.
# Helper function for clusterSystems
# probably not useful elsewhere

# The distance metric is determined by the distance function that is passed as a parameter.
# The main use of this function is as a helper function for \code{\link{clusterSystems}}.

# Parameters:
# dist_fn       Function to be used to calculate individual pairwise distances.
# genes         A vector of HGNC symbols for which to calculate the similarity of all possible pairwise combinaitons.
# data_source   An object that will be passed to the similarity function in order to calculate similarity of the 2 genes.

# Return: Object of class distance

# Author: Rachel Silverstein

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
