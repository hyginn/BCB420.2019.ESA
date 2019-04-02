# get_best_matches.R

# Purpose: Helper function for clusterSystems
# Not useful elsewhere

get_best_matches <- function (systems, clusts) {
  # For each system, determine which cluster best represents it using the Jaccard index
  best_matches <- character(length = length(systems))
  jaccard_indexes <- numeric(length = length(systems))
  names(best_matches) <- names(systems)
  names(jaccard_indexes) <- names(systems)

  for (i in seq_along(systems)) {
    system <- unlist(systems[[i]])
    curr_max <- 0
    curr_best_match <- 1
    for (j in seq_along(systems)) {
      cluster <- clusts[clusts == j]
      genes <- names(cluster)
      union <- length(union(genes, system))
      intersect <- length(intersect(genes, system))
      jaccardIndex <- intersect / union
      if (jaccardIndex > curr_max) {
        curr_max <- jaccardIndex
        curr_best_match <- j
      }
    }
    best_matches[i] <- curr_best_match
    jaccard_indexes[i] <- curr_max
  }
  return(list(best_matches, jaccard_indexes))
}
