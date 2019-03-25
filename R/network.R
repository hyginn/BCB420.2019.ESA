#
# This file produces gene network graph for genes from the system 'PHALY', 
# it shows the relationship between certain genes and the other genes in a system.
#
#get the STRING edges object
StringURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                    "STRINGedges-2019-03-14.RData")
load(url(StringURL)) # loads STRING edges object
load(url(StringURL)) # loads STRING edges object

if(!require("igraph")){
  install.packages("igraph")
}
library("igraph")

# We can select the gene we want to study from the system 'PHALY',
# first, show all the gene in the system

fetchComponents <- function(sys) {
  # returns a fixed set of symbols.
  # Function stub for development purposes only.
  if (sys == "PHALY") {
    s <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2", 
           "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7", 
           "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP", 
           "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM", 
           "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B", 
           "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN", 
           "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21", 
           "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN", 
           "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6", 
           "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1", 
           "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7", 
           "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", 
           "VPS41", "VTI1B", "YKT6")
  } else {
    s <- ""
  }
  return(s)
}

# 1.query a gene network
gene = "VPS41" #randomly select a gene
# find the genes that is attached to the VPS41 through edges 
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
g = graph(edges=networkdata, directed=F)
# plot the graph
plot(g,edge.arrow.size=.2, vertex.color="red", vertex.size=2, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=0.6, vertex.label.dist=1, edge.curved=0.1)
#
# Next, produce graph with multiple genes
#2.query multi genes network
genes = c("ARF5","SPTBN2","KTF22")
# find all the genes that are attached to the selected genes
geneNetwork2 = STRINGedges[STRINGedges$a %in% genes | STRINGedges$b %in% genes,]
# find the unique genes attach to each selected gene
nodes = unique(c(geneNetwork2$a,geneNetwork2$b))
# create igraph graphs from geneNetwork2
net <- graph_from_data_frame(d=geneNetwork2,vertices = nodes,directed = F)
# plot the graph, set the color and length
plot(net, edge.arrow.size=.2, vertex.color=c("red"), vertex.size=2, 
     vertex.frame.color="gray", vertex.label.color="black", 
     vertex.label.cex=0.6, vertex.label.dist=1, edge.curved=0)

# [END]