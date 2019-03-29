#=============================================TOC====================================================
#======================================Section 1: Dataset I choose from GENE POLR1A profile==========
#======================================Section 2: Dataset Provided by professor======================


#============Section 1 is the dataset I from POLR1A profile==========================================
# This partial code was obtained from BCB410 Course.
source("https://bioconductor.org/biocLite.R")
if (!require(Biobase, quietly=TRUE)) {
  biocLite("Biobase")
  library(Biobase)
}
if (!require(GEOquery, quietly=TRUE)) {
  biocLite("GEOquery")
  library(GEOquery)
}
# Load series and platform data from GEO
raw_dataset_1 <- getGEO("GSE34512", GSEMatrix =TRUE, getGPL=FALSE)
if (length(raw_dataset_1) > 1) {
  index_of_raw_dataset_1 <- grep("GPL570", attr(raw_dataset_1, "names"))
} else {
  index_of_raw_dataset_1 <- 1
}
raw_dataset_1 <- raw_dataset_1[[index_of_raw_dataset_1]]
# let us save file into directory.
GSE34512 <- raw_dataset_1
save(GSE34512, file="./data/GSE34512.RData")
load(file="./data/GSE34512.RData")
raw_dataset_1  <-GSE34512
raw_dataset_1
# Access contents via methods:
featureNames(raw_dataset_1 )[1:10]
sampleNames(raw_dataset_1 )
exprs(raw_dataset_1)
source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
if (!require(Biobase, quietly=TRUE)) {
  biocLite("Biobase")
  library(Biobase)
}
if (!require(GEOquery, quietly=TRUE)) {
  biocLite("GEOquery")
  library(GEOquery)
}
# install package outliers to detect any outliers of dataframe
install.packages("outliers")
library(outliers)
# install package class for clustering purpose
install.packages("class")
library(class)
number_of_genes <- nrow(raw_dataset_1)
number_of_genes
number_of_samples <- ncol(raw_dataset_1)
number_of_samples
# assign all of features and columns into an object;
expression_dataset_1 <- raw_dataset_1[1:10,1:number_of_samples]
expression_dataset_1
# assign above object into a matrix using "exprs" function
expression_dataset_1_Matrix <- matrix(exprs(expression_dataset_1), nrow=10,ncol=number_of_samples)
expression_dataset_1_Matrix
ncol(expression_dataset_1_Matrix)
nrow(expression_dataset_1_Matrix)
# omit the rows that has na value;
na.omit(expression_dataset_1_Matrix)
ncol(expression_dataset_1_Matrix)
nrow(expression_dataset_1_Matrix)
#lets transpose the matrix
expression_dataset_1_Matrix_transpose <- t(expression_dataset_1_Matrix)
expression_dataset_1_Matrix_transpose
# let's remove the outliers
rm.outlier(expression_dataset_1_Matrix_transpose)
nrow1 <- nrow(expression_dataset_1_Matrix_transpose)
nrow1
ncol1 <- ncol(expression_dataset_1_Matrix_transpose)
ncol1
# create correlation matrix
expression_dataset_1_Matrix_transpose_cor <- cor(expression_dataset_1_Matrix_transpose)
expression_dataset_1_Matrix_transpose_cor
threshold_value <- sequence(0.1:1)
# Let's create adjancency matrix by choosing threhold value to 0.3.
expression_dataset_1_Matrix_transpose_cor_adjancency_0.3 <- expression_dataset_1_Matrix_transpose_cor

for (c in 1:nrow(expression_dataset_1_Matrix_transpose_cor))
{
  for (d in 1:ncol(expression_dataset_1_Matrix_transpose_cor))
  {
    if(expression_dataset_1_Matrix_transpose_cor[c,d] >= 0.3)
    {
      expression_dataset_1_Matrix_transpose_cor_adjancency_0.3[c,d] <-1
    }
    else {
      expression_dataset_1_Matrix_transpose_cor_adjancency_0.3[c,d] <-0
    }
  }

}
for (c in 1:nrow(expression_dataset_1_Matrix_transpose_cor))
{
  for (d in 1:ncol(expression_dataset_1_Matrix_transpose_cor))
  {
    if(expression_dataset_1_Matrix_transpose_cor[c,d] >= 0.3)
    {
      expression_dataset_1_Matrix_transpose_cor_adjancency_0.3[c,d] <-1
    }

    else if(expression_dataset_1_Matrix_transpose_cor[c,d] <= -0.3)
    {
      expression_dataset_1_Matrix_transpose_cor_adjancency_0.3[c,d] <- -1
    }
  }

}
expression_dataset_1_Matrix_transpose_cor_adjancency_0.3
install.packages("network")
library("network")
coexpressedgenes <- network(expression_dataset_1_Matrix_transpose_cor_adjancency_0.3, directed=FALSE)
# plot the the network of these graph
gene_name <- featureNames(raw_dataset_1 )[1:10]
plot(coexpressedgenes, label=gene_name)
#=======================Igraph ploting================================================================
install.packages("igraph")
library("igraph")
# let's create an igrpah object from previous adjancency matrix.
coexpressed_genes<-graph_from_adjacency_matrix(expression_dataset_1_Matrix_transpose_cor_adjancency_0.3
                                               , mode="undirected",weighted=TRUE)
igraph_1 <- simplify(coexpressed_genes)
V(igraph_1)$name<-gene_name
plot.igraph(igraph_1,vertex.color="red", vertex.size = 10,directed=FALSE)
# This is non-interactive 2D plot meaning you can't rotate the graph. we can use "tkplot" function
#for interactive 2D plot.
#If you are running OS X, you need to install XQuartz which can be downlaod at https://www.xquartz.org,
tkplot(igraph_1) # this is interactive 2D plotting.
#if we want to see 3D ploting, we need to install "rgl" package
intall.packages("rgl")
library(rgl)
tkplot(igraph_1)
rglplot(igraph_1) #this will give you 3D plot
#if we want to find out co-expressed genes, we can just enter gene name and find out the co-expressed
#genes.
gene <- c("117_at")
neighbours <- neighbors(igraph_1,"117_at",mode="all")
neighbours
#   + 2/10 vertices, named, from b4e6afc:
#  [1] 1255_g_at 1431_at
subgraph_1255_g_at <-c("117_at","1255_g_at", "1431_at")
# Let's plot coexpressed genes;
plot.igraph(induced_subgraph(igraph_1,subgraph_1255_g_at))
# let's find out either the graph is star-like or not by calculating the centrality
(centrality <- centralization.betweenness(igraph_1, directed=TRUE, nobigint =TRUE,normalized = TRUE))
#The centralization score is only 0.3796296(rougly 38%) which is not really centralzied.
#cliques in graph means every vertices in graph are adjacent and means complete graph. In gene expression,
# it means cluster and these protein interact with each other closely.
(largest_cliques(igraph_1))
#[[1]]
#+ 6/10 vertices, named, from b4e6afc:
#  [1] 121_at    1007_s_at 1053_at   1294_at   1320_at   1405_i_at
#it returns the vertices with the subgraph with highest vertices. This could mean these genes have
#similar functoionality.
(max_cliques(igraph_1))
cluster_edge_betweenness(igraph_1,weights = NULL)# it returns a community which is similar to large clique
#IGRAPH clustering edge betweenness, groups: 3, mod: 0.25
#+ groups:
#  $`1`
#[1] "1007_s_at" "1053_at"   "121_at"    "1294_at"   "1320_at"   "1405_i_at"
#$`2`
#[1] "117_at"    "1255_g_at" "1431_at"
#$`3`
#[1] "1316_at"


#==================================Section 2: Dataset Provided by Professor==========================
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "GEO-QN-profile-2019-03-24.rds")
myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data

str(myQNXP)


sum( ! is.na(myQNXP))

mean(  myQNXP, na.rm = TRUE)


median(myQNXP, na.rm = TRUE)

hist(log10(myQNXP), breaks = 100)

sum( ! is.na(rowMeans(myQNXP, na.rm = TRUE)))

sum( ! is.na(rowMeans(myQNXP, na.rm = FALSE)))

colnames(myQNXP)

myExp <- gsub("^([^\\.]+).+", "\\1", colnames(myQNXP))

length(unique(myExp))

myQNXP["VAMP8", ]


# plot expression values
plot(myQNXP["VAMP8", ], ylab = "quantile normalized expression (AU)")

# add lines that separate the experiments
abline(v = cumsum(rle(myExp)$lengths) + 0.5, lwd = 0.5, col = "#0000CC44")

```

Gene / Gene correlations
```R
corGenes <- function(A, B, prf) {
  # Calculate pearson correlation between gene expression
  # profiles A and B in prf identified by the gene symbol.
  # A and B can be either gene symbol or index.

  r <- cor(prf[A, ], prf[B, ], use = "pairwise.complete.obs")
  return(r)
}

corGenes("MRPL18", "RPF2", myQNXP)
corGenes(100, 200, myQNXP)

plotCorGenes <- function(A, B, prf) {
  # Plot correlation between gene expression
  # profiles A and B in prf identified by the gene symbol.
  # A and B can be either gene symbol or index.

  xMin <- min(c(0, prf[A, ], prf[B, ]), na.rm = TRUE)
  xMax <- max(c(   prf[A, ], prf[B, ]), na.rm = TRUE)

  plot(prf[A, ], prf[B, ],
       xlim = c(xMin, xMax),
       ylim = c(xMin, xMax),
       main = sprintf("%s vs. %s (r = %5.2f)", A, B, corGenes(A, B, prf)),
       xlab = (sprintf("%s expresion (AU)", A)),
       ylab = (sprintf("%s expresion (AU)", B)),
       asp = 1.0)
  abline(lm(prf[B, ] ~ prf[A, ]), col = "#AA0000", lwd = 0.67)

  return(invisible(NULL))
}

plotCorGenes(A = "SLC9A4", B = "NLGN2", prf = myQNXP)






