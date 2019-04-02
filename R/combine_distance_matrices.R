# combine_distance_matrices.R

#' Combine Distance Matrices
#'
#' \code{combine_distance_matrices} This is a helper function for \code{\link{clusterSystems}} and is unlikely to be useful elsewhere.
#'
#' @param mode Method by which to combine distance matrices
#' @param distanceMatrices List of objects of class distance.
#'
#' @return single object of class distance.
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#'
combine_distance_matrices <- function(mode, distanceMatrices) {
  if (length(distanceMatrices) > 1) {

    if (mode == 'sum') {
      for (i in seq_along(distanceMatrices)) {
        if (i == 1) {
          combinedMatrix <- distanceMatrices[[i]]
        } else {
          combinedMatrix <- combinedMatrix + distanceMatrices[[i]]
        }
      }
    } else if (mode == 'product') {
      for (i in seq_along(distanceMatrices)) {
        if (i == 1) {
          combinedMatrix <- distanceMatrices[[i]]
        } else {
          combinedMatrix <- combinedMatrix * distanceMatrices[[i]]
        }
      }
    } else if (mode == 'minimum') {
      for (i in seq_along(distanceMatrices)) {
        if (i == 1) {
          combinedMatrix <- distanceMatrices[[i]]
        } else {
          combinedMatrix <- pmin(combinedMatrix,  distanceMatrices[[i]])
        }
      }
    } else if (mode == 'maximum') {
      for (i in seq_along(distanceMatrices)) {
        if (i == 1) {
          combinedMatrix <- distanceMatrices[[i]]
        } else {
          combinedMatrix <- pmax(combinedMatrix, distanceMatrices[[i]])
        }
      }
    }

  } else {
    combinedMatrix <- distanceMatrices[[1]]
  }
  return(combinedMatrix)
}
