# make_venn_diagram.R

#' Make Venn Diagrams
#'
#' \code{make_venn_diagram} is a helper function for \code{\link{clusterSystems}} and is unlikely to be useful elsewhere.
#'
#' @param system_genes vector of genes in the system
#' @param sys_name name of the system
#' @param cluster_genes vector of genes in the cluster
#' @param clust_num Numeric identifier for the cluster
#' @param jaccard The jaccard index of the set overlap of cluster and system
#' @param p P value of the Jaccar index
#'
#' @return NULL
#'
#' @import VennDiagram
#'
#' @author \href{https://orcid.org/0000-0001-5724-2252}{Rachel Silverstein} (aut)
#'
#'



make_venn_diagram <- function(system_genes, sys_name, cluster_genes, clust_num, jaccard, p) {

  clust_name <- paste(c("Cluster ", clust_num), collapse = "")


  area1 <- length(system_genes)
  area2 <- length(cluster_genes)
  int <- length(intersect(system_genes, cluster_genes))
  category <- c(sys_name, clust_name)

  grid::grid.newpage()
  venn <- VennDiagram::draw.pairwise.venn(area1 = area1,
                                            area2 = area2,
                                            cross.area = int,
                                            category = category,
                                            euler.d = T,
                                            scaled = T,
                                            fill = c('red', 'blue'),
                                            alpha = c(0.5, 0.5),
                                            cat.default.pos = 'outer')

  return(NULL)
}

