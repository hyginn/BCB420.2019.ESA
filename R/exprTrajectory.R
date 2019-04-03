# exprTrajectory.R

#' Title.
#'
#' \code{plotExprTrajectory} Plots an iGraph of the p-Creode with coloring of the gene expression
#'
#' Details.
#' @section Dependency: the pcreode package, which can be installed with:
#' \code{sudo pip install pcreode}
#'
#' @param overlayName The specific Gene/Protein in the .fcs file you would like to color on the trajectories
#' @return
#'
#' @author Yoonsik Park (au)
#'
#' @examples
#' # plots the CAR2 gene expression for the current expression trajectory
#' plotExprTrajectory("CAR2")
#'
#' @export

plotExprTrajectory <- function (overlayName) {
  system(paste("/usr/local/bin/python ./dev/plotExprTrajectory.py ./data/temp/temp.csv ", overlayName, sep="", collapse = NULL))

  m <- as.matrix(read.table("./data/temp/weights.txt"))
  dens <- as.matrix(read.table("./data/temp/density.txt"))
  nodes <- as.matrix(read.table("./data/temp/nodes.txt"))
  colors <- as.matrix(read.table("./data/temp/colors.txt"))

  ig <- igraph::graph.adjacency(m, mode="undirected", weighted=TRUE)
  layout <- igraph::layout.kamada.kawai(ig, maxiter=10000)
  graphCol = c()
  for (i in 1:nrow(colors[,1:3])) {
    graphCol = c(graphCol,colorspace::hex(colorspace::RGB(t(as.matrix(colors[i,3:1])))))
  }
  titleString = paste("p-Creode of Expression Data colored by ", overlayName, sep="", collapse = NULL)
  plot(ig, vertex.size=dens, edge.width=2,vertex.color=graphCol, vertex.label=NA, layout=layout, main=titleString)
}

#' Title.
#'
#' \code{createExprTrajectory} Creates the p-Creode Graphs from flow cytometry
#'
#' Details.
#' @section Dependency: the pcreode package, which can be installed with:
#' \code{sudo pip install pcreode}
#'
#' @param filename The filename of the .fcs flow cytometry file.
#' @return
#'
#' @author Yoonsik Park (au)
#'
#' @examples
#' # Loads the Myeloid Differentiation Flow Cyt. Dataset, and stores the trajectory graphs
#' createExprTrajectory("./inst/extdata/rna_seq_myeloid.fcs")
#'
#' @export

createExprTrajectory <- function (filename) {
  pre_data <- flowCore::read.FCS(filename)
  write.csv(data.frame(pre_data@exprs), "./data/temp/temp.csv", row.names=F)
  system("/usr/local/bin/python ./dev/createExprTrajectory.py ./data/temp/temp.csv ")
  }

# [END]
