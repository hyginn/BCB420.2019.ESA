# make_venn_diagram.R

# Purpose:  a helper function for clusterSystems and is unlikely to be useful elsewhere.

# Parameters:
# system_genes      vector of genes in the system
# sys_name          name of the system
# cluster_genes     vector of genes in the cluster
# clust_num         Numeric identifier for the cluster
#
# Value: NULL

# Author: Rachel Silverstein
#' @import grid

make_venn_diagram <- function(system_genes, sys_name, cluster_genes, clust_num) {

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

