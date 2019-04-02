# ClusterSystems.R

#' Cluster Systems
#'
#' \code{ClusterSystems} inputs a list of biological systems and re-clusters
#' the genes on an input variable of interest. Set overlap
#' between the output clusters and clusters defined by input system
#' membership is measured.
#'
#' The input systems are clustered according to the specified variables(s)
#' of interest using PAM (partition around medoids) clustering. The
#' similarity between input and output sets is measured using the Jaccard
#' index of set overlap. P values for the observed Jaccard Indexes are
#' calculated by measuring the Jaccard index of 1000 clusters randomly
#' sampled from the input genes.
#'
#' @section Distance metrics :
#' There are several different built in methods to compute the distance
#' between a pair of genes.
#'
#' \code{"expression_profile"} uses GEO data to compute the similarity
#' between 2 genes as the absolute value of the Spearmann correlation
#' between their expression profiles. Distance is then taken as 1 - similarity.
#'
#' \code{"transcription_factor"} GTRD data to compute the number of
#' shared transcription that bind upstream of the 2 genes divided by the
#' total number of transcription factors.
#' In effect, the Jaccard index of the sets of transcription factors from the
#' 2 genes. Distance is then calculated as 1 - similarity.
#'
#' \code{"network_jaccard"} uses STRING network data to calculate the Jaccard
#' similarity between the immediate neighbours of the 2 genes. Distance is
#' then calculated as 1 - similarity.
#'
#' \code{"network_distance"} uses STRING network data and calculates the
#' shortest path between the 2 genes.
#'
#' @section More than one distance metric: If more than one distance data
#' type is provided, the various distance metrics will be combined into one
#' distance matrix according to one of the following methods:
#'
#' \code{'sum'} indicates the distances between a pair of genes will be
#' defined as the sum of the distances according to each metric.
#'
#' \code{'product'} indicates that the distance between genes is defined as
#' the product of their pairwise distances. This method more strongly
#' penalizes genes that are distant by more than one metric.
#'
#' \code{'maximum'} or \code{'minimum'} indicate that the distance between
#' 2 genes should be the maximum or minimum of the distances measured by each
#' metric respectively.
#'
#' @section Custom distance metrics: In addition to several distance metrics
#' which are built in to the function, the user has the option of defining a
#' new distance metric to measure pairwise similarity between genes. If this
#' is the case, an appropriate source of data (mapping from HGNC symbols to
#' data of interest) that can be used as in input for each custom distance
#' function.
#'
#' @param systems List of input systems that should be
#' re-clustered. Each element is formatted according to the output of
#' \code{SyDBgetSysSymbols} (a named list of HGNC symbols).
#' @param distances Character vector indicating choice of data to be used
#' to measure distance between input genes. Vector containing one or more of
#'  \code{c("expression_profile", "transcription_factor", "network_jaccard",
#'  "network_distance")} or NULL if no built in distance metric is to be used
#'  (in this case a custom distance metric must be used instead.) Default is NULL.
#'
#' @param customDistanceFn Optional list of functions. Any custom functions
#' that are to be used to calculate pairwise distances between genes. Default
#' is NULL. If this parameter is not NULL, then \code{dataSources} parameter
#' must also be provided and be a list of the same length.
#' @param dataSources Optional list of input data appropriate to be used in
#' conjunction with \code{customDistanceFn}. Elements may have any type
#' that is required as input by the custom distance function, but must be
#' contained within a list. The entries should
#' be in the same order as \code{customDistanceFn} such that the ith element
#' of \code{customDistanceFn} uses the data provided in the ith element of
#' \code{dataSources} as input. Default is NULL.
#' @param combineMatrices Length 1 character vector indicating the method by which
#' to combine distance matrices if more than one choice of distance data is provided.
#' If only one type of distance data is provided, then this parameter is
#' ignored. One of \code{c('sum', 'product', 'minimum', 'maximum')}.
#' @param plotVennDiagrams Logical flag; indicated whether or not Venn
#' diagrams representing the set overlap of the output should be printed when
#' the function is called. Default is TRUE.
#' @param k Optional integer. The number of clusters. Default is \code{length(systems)}
#' if no value provided.
#'
#' @return A named list of length 4. Elements of the list include:
#'
#' \code{$Clusters} A named vector of integers from 1 to k. Names of the
#' elements represent the genes belonging to that cluster.
#'
#' \code{$Best_matches} A named vector of length k representing the cluster
#' with the best set overlap for each input system.
#'
#' \code{$Jaccard_indexes} The Jaccard indexes measuring the similarity
#' between each system and the cluster which is its best match.
#'
#' \code{$P_values} P values representing the probability of having a
#' Jaccard index greater than those observed by choosing a random cluster
#' of the same size from the set of input genes.
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#' @seealso For examples of distance functions that can be used to create
#' custom distance matrices, see \code{\link{expr_dist}},
#' \code{\link{jaccard_dist}}, and \code{\link{tf_dist}}.
#'
#' @examples
#'
#' myDB <- fetchData("SysDB")
#' rootSysIDs <- SyDBgetRootSysIDs(myDB)
#' sys_names <- names(rootSysIDs)
#' systems <- SyDBgetSysSymbols(myDB, sys_names)
#'
#' # Cluster all of the systems in the database according to
#' # the sum of transcription factor distance and Jaccard network distance:
#' clusterSystems(systems,
#'                distances = c("transcription_factor", "network_jaccard"),
#'                combineMatrices = 'sum')
#'
#' # Cluster all of the systems using expression profile distance
#' # but format it like a custom distance function
#' GEO <- fetchData("GEOprofiles")
#' clusterSystems(systems,
#'                distances = NULL,
#'                customDistanceFn = list(expr_dist),
#'                dataSources = list(GEO))
#'
#'
#'
#' @import cluster
#' @import VennDiagram
#' @import igraph
#'
#' @export


clusterSystems <- function(systems,
                           distances = NULL,
                           customDistanceFn = NULL,
                           dataSources = NULL,
                           combineMatrices,
                           plotVennDiagrams = TRUE,
                           k) {

  # check to make sure at least one distance metric is provided
  if (is.null(distances) & is.null(customDistanceFn)) {
    stop("A value for at least one of <distances> or <customDistanceFn> must be provided.")
  }

  # make a single vector containing all of the genes from all of the systems together
  all_genes <- unique(unlist(systems))

  distanceMatrices <- make_distance_matrices_list(distances, all_genes)

  # compute any distance matrices using custom similarity functions and data sources
  if (!is.null(customDistanceFn)) {
    for (i in seq_along(customDistanceFn)) {
      data <- dataSources[[i]]
      fn <- customDistanceFn[[i]]
      matrix <- make_matrix(all_genes, fn, data)
      matrix <- list(matrix)
      names(matrix) <- paste(c("custom_function_", i), collapse = "")
      distanceMatrices <- append(distanceMatrices, matrix)
    }
  }

  # Combine the distance matrices if there is more than one
  combinedMatrix <- combine_distance_matrices(mode = combineMatrices,
                                              distanceMatrices)


  # CLUSTER TIME :-)

  # if k missing use the number of systems as the number of clusters for PAM clustering
  if (missing(k)) { k <- length(systems) }

  clust_obj <- cluster::pam(combinedMatrix, k = k, diss = TRUE)

  clusters <- clust_obj$clustering

  # get the best cluster - system pairings
  matches <- get_best_matches(systems, clusters)
  best_matches <- matches[[1]]
  jaccard_indexes <- matches[[2]]

  # Calculate P values of getting Jaccard index this high or higher for each
  # system - cluster pair
  p_values <- jaccard_index_p_values(systems,
                                     jaccard_indexes,
                                     best_matches,
                                     clusters,
                                     trials = 1000)

  result <- list("Clusters" = clusters,
                 "Best_matches" = best_matches,
                 "Jaccard_indexes" = jaccard_indexes,
                 "P_values" = p_values)

  # display the results as venn diagrams
  if (plotVennDiagrams == TRUE) {
    for (i in seq_along(systems)) {

      clust_num <- best_matches[i]
      sys_name <- names(systems)[i]
      cluster_genes <- names(clusters[clusters == clust_num])
      system_genes <- unlist(systems[[i]])
      jaccard <- jaccard_indexes[i]
      p <- p_values[i]
      make_venn_diagram(system_genes, sys_name, cluster_genes, clust_num)
    }
  }

  return(result)
}

