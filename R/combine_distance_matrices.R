# combine_distance_matrices.R

# This is a helper function for clusterSystems and is unlikely to be useful elsewhere.

# Parameters :
# mode                Method by which to combine distance matrices
# distanceMatrices    List of objects of class distance.
#
# Value:
# single object of class distance.
#
# Author: Rachel Silverstein

combine_distance_matrices <- function(mode, distanceMatrices) {
  if (length(distanceMatrices) > 1) {

    for (i in seq_along(distanceMatrices)) {
      if (i == 1) {
        combinedMatrix <- distanceMatrices[[i]]
      } else {
        if (mode == 'sum') {
          combinedMatrix <- combinedMatrix + distanceMatrices[[i]]
        } else if (mode == 'product') {
          combinedMatrix <- combinedMatrix * distanceMatrices[[i]]
        } else if (mode == 'maximum') {
          combinedMatrix <- pmax(combinedMatrix, distanceMatrices[[i]])
        } else if (mode == 'minimum') {
          combinedMatrix <- pmin(combinedMatrix,  distanceMatrices[[i]])
        }
      }
    }

    # if (mode == 'sum') {
    #   for (i in seq_along(distanceMatrices)) {
    #     if (i == 1) {
    #       combinedMatrix <- distanceMatrices[[i]]
    #     } else {
    #       combinedMatrix <- combinedMatrix + distanceMatrices[[i]]
    #     }
    #   }
    # } else if (mode == 'product') {
    #   for (i in seq_along(distanceMatrices)) {
    #     if (i == 1) {
    #       combinedMatrix <- distanceMatrices[[i]]
    #     } else {
    #       combinedMatrix <- combinedMatrix * distanceMatrices[[i]]
    #     }
    #   }
    # } else if (mode == 'minimum') {
    #   for (i in seq_along(distanceMatrices)) {
    #     if (i == 1) {
    #       combinedMatrix <- distanceMatrices[[i]]
    #     } else {
    #       combinedMatrix <- pmin(combinedMatrix,  distanceMatrices[[i]])
    #     }
    #   }
    # } else if (mode == 'maximum') {
    #   for (i in seq_along(distanceMatrices)) {
    #     if (i == 1) {
    #       combinedMatrix <- distanceMatrices[[i]]
    #     } else {
    #       combinedMatrix <- pmax(combinedMatrix, distanceMatrices[[i]])
    #     }
    #   }
    # }

  } else {
    combinedMatrix <- distanceMatrices[[1]]
  }
  return(combinedMatrix)
}
