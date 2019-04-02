# network.R

#' \code{produceGeneNetwork} produces gene network graph for genes from the given system of genes 'sys'.
#'
#' \code{produceGeneNetwork} outputs a graph on a network graph corresponding to a selected gene
#' in a system, it shows the correlation between one gene to the rest in a system.
#'
#' @param sys (character)       5-character system code of genes
#' @param gene (character)      genes code

#' @return a network graph
#' @import igraph
#' @examples
#' produceGeneNetwork("PHALY", "VPS41")

#' @export

produceGeneNetwork <- function(sys, gene="VPS41") {
  StringURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                      "STRINGedges-2019-03-14.RData")
  load(url(StringURL)) # loads STRING edges object
  load(url(StringURL)) # loads STRING edges object

  # fetch data
  myDB <- fetchData("SysDB")
  geneSet <- SyDBgetSysSymbols(myDB, sys)

  StringURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                      "STRINGedges-2019-03-14.RData")

  load(url(StringURL)) # loads STRING edges object
  load(url(StringURL)) # loads STRING edges object

  # 1.query a gene network
  # gene = "VPS41" #randomly select a gene
  # find the genes that is attached to the VPS41 by edges
  geneNetwork = STRINGedges[STRINGedges$a == gene|STRINGedges$b==gene ,]

  # get the number of genes in the list
  rowcount = dim(geneNetwork)[1]
  networkdata = c()
  # append all the genes to the list
  for(i in 1:rowcount){
    networkdata = append(networkdata,geneNetwork[i,]$a)
    networkdata = append(networkdata,geneNetwork[i,]$b)
  }
  # add edges between genes
  g = igraph::graph(edges=networkdata, directed=F)
  # plot the graph
  graphics::plot(g,edge.arrow.size=.2, vertex.color="red", vertex.size=2,
       vertex.frame.color="gray", vertex.label.color="black",
       vertex.label.cex=0.6, vertex.label.dist=1, edge.curved=0.1)
}
#' \code{produceGenesNetwork} produces genes network graph for genes from the given system of genes 'sys'.
#'
#' \code{produceGenesNetwork} outputs a graph on a network graph corresponding to a selected gene
#' in a system, it shows the correlation between one gene to the rest in a system.
#'
#' @param sys (character)         5-character system code of genes
#' @param genes (character list)  genes code

#' @return a network graph
#' @import igraph
#' @examples
#' produceGenesNetwork("PHALY", c("ARF5","SPTBN2","KTF22"))

#' @export
produceGenesNetwork <- function(sys, genes=c("ARF5","SPTBN2","KTF22")) {
  StringURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                      "STRINGedges-2019-03-14.RData")
  load(url(StringURL)) # loads STRING edges object
  load(url(StringURL)) # loads STRING edges object
  #fetch data
  myDB <- fetchData("SysDB")
  geneSet <- SyDBgetSysSymbols(myDB, sys)

  StringURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                      "STRINGedges-2019-03-14.RData")
  load(url(StringURL)) # loads STRING edges object
  load(url(StringURL)) # loads STRING edges object
  # find all the genes that are attached to the selected genes
  geneNetwork = STRINGedges[STRINGedges$a %in% genes | STRINGedges$b %in% genes,]
  # find the unique genes attach to each selected gene
  nodes = unique(c(geneNetwork$a,geneNetwork$b))
  # create igraph graphs from geneNetwork2
  net <- igraph::graph_from_data_frame(d=geneNetwork,vertices = nodes,directed = F)
  # plot the graph, set the color and length
  graphics::plot(net, edge.arrow.size=.2, vertex.color=c("red"), vertex.size=2,
       vertex.frame.color="gray", vertex.label.color="black",
       vertex.label.cex=0.6, vertex.label.dist=1, edge.curved=0)
}
# [END]

