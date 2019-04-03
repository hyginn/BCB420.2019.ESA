# `BCB420.2019.ESA`

#### (**E**xploratory **S**ystems **A**nalysis **Tools for BCB420-2019**)
<!-- [![DOI](https://zenodo.org/badge/157482801.svg)](https://zenodo.org/badge/latestdoi/157482801) -->

&nbsp;

###### [Boris Steipe](https://orcid.org/0000-0002-1134-6758),
###### Department of Biochemistry and Department of Molecular Genetics,
###### University of Toronto
###### Canada
###### &lt; boris.steipe@utoronto.ca &gt;

&nbsp;
###### [Nada Elnour](https://orcid.org/0000-0001-6165-1542); nada.elnour@mail.utornoto.ca ######

----

This README describes this package.
<!-- The associated Vignette can be previewed [here](http://htmlpreview.github.io/?https://github.com/hyginn/rptPlus/blob/master/doc/rptPlusVignette.html). The package can be installed from GitHub with `devtools::install_github("hyginn/rptPlus", build_opts = c("--no-resave-data", "--no-manual"))`. -->

----

**If any of this information is ambiguous, inaccurate, outdated, or incomplete,
please check the [most recent version](https://github.com/hyginn/BCB420.2019.ESA) of the
package on GitHub and 
[file an issue](https://github.com/hyginn/BCB420.2019.ESA/issues).**

----

<!-- TOCbelow -->
1. About this package:<br/>
2. Data ...<br/>
3. Functions ... <br />
4. Notes<br/>
5. References and Further reading<br/>
6. Acknowledgements<br/>
<!-- TOCabove -->

----

## 1 About this package:

This package is a joint development platform for student designed tools for Exploratory Systems Analysis, part of the University of Toronto course BCB420H1S (Computational Systems Biology) in the 2018/2019 academic year.

&nbsp;


----

## 2 Data ...

Suporting resources include curated systems data and other data resources:

#### 2.1 The HGNC symbol reference

Load the `HGNC` object in the following way:

```R
myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame

str(HGNC)
# 'data.frame':	27087 obs. of  14 variables:
#  $ sym      : chr  "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#  $ name     : chr  "alpha-1-B glycoprotein" "A1BG antisense RNA 1"  ...
#  $ UniProtID: chr  "P04217" NA "Q9NQ94" "P01023" ...
#  $ RefSeqID : chr  "NM_130786" "NR_015380" "NM_001198818" "NM_000014" ...
#  $ EnsID    : chr  "ENSG00000121410" "ENSG00000268895" "ENSG00000148584" ...
#  $ UCSCID   : chr  "uc002qsd.5" "uc002qse.3" "uc057tgv.1" "uc001qvk.2" ...
#  $ GeneID   : num  1 503538 29974 2 144571 ...
#  $ OMIMID   : chr  "138670" NA "618199" "103950" ...
#  $ acc      : chr  NA "BC040926" "AF271790" "BX647329, X68728, M11313" ...
#  $ chr      : chr  "19q13.43" "19q13.43" "10q11.23" "12p13.31" ...
#  $ type     : chr  "protein" "lncRNA" "protein" "protein" ...
#  $ prev     : chr  NA "NCRNA00181, A1BGAS, A1BG-AS" NA NA ...
#  $ synonym  : chr  NA "FLJ23569" "ACF, ASP, ACF64, ACF65, APOBEC1CF"  ...
#  $ RefSeqOld: chr  "NM_130786" "NR_015380" "NM_014576" "NM_000014" ...

```

&nbsp;

#### 2.2 Transcription factors from GTRD:

`geneList` list contains transcription factors that have been found to bind to the upstream regulatory regions of genes in ChIP-seq experiments; the data is compiled by the GTRD database.
Load `geneList` in the following way:

```R
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "geneList-2019-03-13.RData")
load(url(myURL))  # loads GTRD geneList object

str(geneList)
# List of 17864
#  $ HES4        : chr [1:82] "AR" "ATF2" "ATF3" "ATF4" ...
#  $ PUSL1       : chr [1:64] "AR" "BHLHE40" "CEBPA" "CEBPB" ...
#  $ ACAP3       : chr [1:69] "BHLHE40" "CEBPA" "CEBPB" "CEBPD" ...
#  $ ATAD3B      : chr [1:110] "AR" "ASCL2" "ATF2" "ATF3" ...
#  $ ATAD3A      : chr [1:128] "ARID4B" "ASCL2" "ATF1" "ATF3" ...
#   [list output truncated]

```

&nbsp;

#### 2.3 STRING edges:

`STRINGedges` contains edges of the STRING database mapped to HGNC symbols.
Load `STRINGedges` in the following way:

```R

myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "STRINGedges-2019-03-14.RData")
load(url(myURL))  # loads STRING edges object

str(STRINGedges)
#  Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	319997 obs. of  3 variables:
#   $ a    : chr  "ARF5" "ARF5" "ARF5" "ARF5" ...
#   $ b    : chr  "SPTBN2" "KIF13B" "KIF21A" "TMED7" ...
#   $ score: num  909 910 910 906 971 915 905 927 914 902 ...
```

&nbsp;


#### 2.4 Expression profiles:

Expression profiles were compiled from 52 microarray experiments downloaded from GEO, and quantile normalized. Details to follow. Here is code to load and use the profiles.

```R

# Load the expression profiles:
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "GEO-QN-profile-2019-03-24.rds")
myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data

str(myQNXP)
#  num [1:27087, 1:52] 29.4 199.9 34.3 947.4 2249.2 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:27087] "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#   ..$ : chr [1:52] "GSE35330.ctrl.4h" "GSE35330.cond.4h" "GSE35330.ctrl.16h" ...


# Some statistics:
sum( ! is.na(myQNXP)) # Number of measurements: 923361
mean(  myQNXP, na.rm = TRUE) # 502.7972
median(myQNXP, na.rm = TRUE) # 127.43
hist(log10(myQNXP), breaks = 100)

# how many HGNC genes have at least one measurement recorded?
sum( ! is.na(rowMeans(myQNXP, na.rm = TRUE))) # 21063

# how many HGNC genes have measurements recorded an all experiments?
sum( ! is.na(rowMeans(myQNXP, na.rm = FALSE))) # 10812

# colnames describe experiments (averaged over replicates)
colnames(myQNXP)

# how many unique experiments
myExp <- gsub("^([^\\.]+).+", "\\1", colnames(myQNXP))
length(unique(myExp)) # 15

#Access one profile by name
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

# Load the expression profiles:
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "GEO-QN-profile-2019-03-24.rds")
myQNXP <- readRDS(url(myURL))  # loads quantile-normalized expression data

str(myQNXP)
#  num [1:27087, 1:52] 29.4 199.9 34.3 947.4 2249.2 ...
#  - attr(*, "dimnames")=List of 2
#   ..$ : chr [1:27087] "A1BG" "A1BG-AS1" "A1CF" "A2M" ...
#   ..$ : chr [1:52] "GSE35330.ctrl.4h" "GSE35330.cond.4h" "GSE35330.ctrl.16h" ...


# Some statistics:
sum( ! is.na(myQNXP)) # Number of measurements: 923361
mean(  myQNXP, na.rm = TRUE) # 502.7972
median(myQNXP, na.rm = TRUE) # 127.43
hist(log10(myQNXP), breaks = 100)

# how many HGNC genes have at least one measurement recorded?
sum( ! is.na(rowMeans(myQNXP, na.rm = TRUE))) # 21063

# how many HGNC genes have measurements recorded an all experiments?
sum( ! is.na(rowMeans(myQNXP, na.rm = FALSE))) # 10812

# colnames describe experiments (averaged over replicates)
colnames(myQNXP)

# how many unique experiments
myExp <- gsub("^([^\\.]+).+", "\\1", colnames(myQNXP))
length(unique(myExp)) # 15

#Access one profile by name
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

```
![](./inst/img/lowPairwiseCorrelation.svg?sanitize=true "uncorrelated genes") &nbsp; ![](./inst/img/highPairwiseCorrelation.svg?sanitize=true "highly correlated genes")<br/>
Example plots for uncorrelated and highly correlated gene expression profiles.

&nbsp;

![](./inst/img/QN-GEOprofileCorrelations.svg?sanitize=true "distribution of expression correlations")<br/>
Distribution of expression correlations between 10,000 randomly chosen gene pairs, compared to the distribution of shuffled expression values for the two genes. The correlations are not symmetric around zero. There are more negative correlations than expected, and more excess negative than positive correlations, but if two genes are positively correlated, the correlation trends to be higher.

&nbsp;

#### 2.5 InterPro domain annotations:

Interpro domain annotations were parsed from the 50 GB annotation source of InterPro v. 73.0 and mapped to HGNC symbols (details to follow).

```R
# load the gene-annotations list:
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "genesIPR.rds")
genesIPR <- readRDS(url(myURL))  # reads and assigns genesIPR list

str(genesIPR, list.len = 5)
# List of 19174
#  $ NUDT4B     : chr [1:3] "IPR000086" "IPR015797" "IPR020084"
#  $ IGLV4-69   : chr [1:5] "IPR003599" "IPR007110" "IPR013106" "IPR013783" ...
#  $ IGLV8-61   : chr [1:5] "IPR003599" "IPR007110" "IPR013106" "IPR013783" ...
#  $ IGLV4-60   : chr [1:5] "IPR003599" "IPR007110" "IPR013106" "IPR013783" ...
#  $ IGLV10-54  : chr [1:4] "IPR007110" "IPR013106" "IPR013783" "IPR036179"
#   [list output truncated]


# load the IPR_domain-locations list:
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "IPRgenes.rds")
IPRgenes <- readRDS(url(myURL))  # reads and assigns IPRgenes list

str(IPRgenes, list.len = 5)
# List of 16178
#  $ IPR000086: chr [1:27] "NUDT4B" "NUDT19" "NUDT21" "TRPM2" ...
#  $ IPR015797: chr [1:28] "NUDT4B" "NUDT19" "NUDT21" "TRPM2" ...
#  $ IPR020084: chr [1:13] "NUDT4B" "NUDT3" "NUDT1" "NUDT2" ...
#  $ IPR003599: chr [1:456] "IGLV4-69" "IGLV8-61" "IGLV4-60" "IGLV7-46" ...
#  $ IPR007110: chr [1:674] "IGLV4-69" "IGLV8-61" "IGLV4-60" "IGLV10-54" ...
#   [list output truncated]

```
![](./inst/img/lowPairwiseCorrelation.svg?sanitize=true "uncorrelated genes") &nbsp; ![](./inst/img/highPairwiseCorrelation.svg?sanitize=true "highly correlated genes")<br/>
Example plots for uncorrelated and highly correlated gene expression profiles.

&nbsp;

![](./inst/img/QN-GEOprofileCorrelations.svg?sanitize=true "distribution of expression correlations")<br/>
Distribution of expression correlations between 10,000 randomly chosen gene pairs, compared to the distribution of shuffled expression values for the two genes. The correlations are not symmetric around zero. There are more negative correlations than expected, and more excess negative than positive correlations, but if two genes are positively correlated, the correlation trends to be higher.


&nbsp;

#### 2.5 InterPro domain annotations:

Interpro domain annotations were parsed from the 50 GB annotation source of InterPro v. 73.0 and mapped to HGNC symbols (details to follow).

```R
# load the gene-annotations list:
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "genesIPR.rds")
genesIPR <- readRDS(url(myURL))  # reads and assigns genesIPR list

str(genesIPR, list.len = 5)
# List of 19174
#  $ NUDT4B     : chr [1:3] "IPR000086" "IPR015797" "IPR020084"
#  $ IGLV4-69   : chr [1:5] "IPR003599" "IPR007110" "IPR013106" "IPR013783" ...
#  $ IGLV8-61   : chr [1:5] "IPR003599" "IPR007110" "IPR013106" "IPR013783" ...
#  $ IGLV4-60   : chr [1:5] "IPR003599" "IPR007110" "IPR013106" "IPR013783" ...
#  $ IGLV10-54  : chr [1:4] "IPR007110" "IPR013106" "IPR013783" "IPR036179"
#   [list output truncated]


# load the IPR_domain-locations list:
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "IPRgenes.rds")
IPRgenes <- readRDS(url(myURL))  # reads and assigns IPRgenes list

str(IPRgenes, list.len = 5)
# List of 16178
#  $ IPR000086: chr [1:27] "NUDT4B" "NUDT19" "NUDT21" "TRPM2" ...
#  $ IPR015797: chr [1:28] "NUDT4B" "NUDT19" "NUDT21" "TRPM2" ...
#  $ IPR020084: chr [1:13] "NUDT4B" "NUDT3" "NUDT1" "NUDT2" ...
#  $ IPR003599: chr [1:456] "IGLV4-69" "IGLV8-61" "IGLV4-60" "IGLV7-46" ...
#  $ IPR007110: chr [1:674] "IGLV4-69" "IGLV8-61" "IGLV4-60" "IGLV10-54" ...
#   [list output truncated]

```


&nbsp;

#### 2.6 Systems:

A systems database (of currently four systems) can be loaded with `fetchData()`. Several utility functions have been added to the package (in `./R/SyDButils.R`). Use `SyDBgetSysSymbols(<database>, <sys>)[[1]]` to access all gene symbols for a single system (or subsystem). Use `unlist(SyDBgetSysSymbols(myDB, SyDBgetRootSysIDs(myDB)), use.names = FALSE)` to access all genes in the database.

```R
myDB <- fetchData("SysDB")
```
#### 2.6 Systems:

A systems database (of currently four systems) can be loaded with `fetchData()`. Several utility functions have been added to the package (in `./R/SyDButils.R`). Use `SyDBgetSysSymbols(<database>, <sys>)[[1]]` to access all gene symbols for a single system (or subsystem). Use `unlist(SyDBgetSysSymbols(myDB, SyDBgetRootSysIDs(myDB)), use.names = FALSE)` to access all genes in the database.

```R
myDB <- fetchData("SysDB")

SyDBgetRootSysIDs(myDB)
#                                        PHALY                                       SLIGR 
#  "give.jams-1d8-648b-1e12-6a9f-65421424affe" "cast.rear-5f1-56ef-1cd2-1ae3-54a5556d59ff" 
#                                        NLRIN                                       HVGCR 
#  "scar.blur-9bc-29cf-31f2-1981-4d92edf4d0e6" "help.mink-f96-e98b-9e12-ab41-8217a3ecb0cd" 


names(SyDBgetRootSysIDs(myDB))
# [1] "PHALY" "SLIGR" "NLRIN" "HVGCR"

SyDBgetSysSymbols(myDB, "HVGCR")  # Note: returns a list.
# $HVGCR
#  [1] "ADCY9"    "ADRB2"    "AKAP7"    "CACNA1C"  "CACNA2D1" "CACNA2D3" "CACNB1"   "CACNB3"  
#  [9] "CACNB4"   "CACNG1"   "CACNG6"   "CALM1"    "CALM2"    "CALM3"    "GNAS"     "PRKACA"  
# [17] "PRKACB"   "PRKAR1A"  "PRKAR1B"  "PRKAR2A"  "PRKAR2B" 

```
&nbsp;

```R

SyDBgetSysSymbols(myDB, "HVGCR")  # Note: returns a list.
# $HVGCR
#  [1] "ADCY9"    "ADRB2"    "AKAP7"    "CACNA1C"  "CACNA2D1" "CACNA2D3" "CACNB1"   "CACNB3"  
#  [9] "CACNB4"   "CACNG1"   "CACNG6"   "CALM1"    "CALM2"    "CALM3"    "GNAS"     "PRKACA"  
# [17] "PRKACB"   "PRKAR1A"  "PRKAR1B"  "PRKAR2A"  "PRKAR2B" 

```
&nbsp;

The old function stub `fetchComponents()` still works but is deprecated:

```R
> fetchComponents("PHALY")
Note: fetchComponents(...) is deprecated. Use SyDBgetSysSymbols(<database>, <system code>) instead.

 [1] "AMBRA1"    "ATG14"     "ATP2A1"    "ATP2A2"    "ATP2A3"    "BECN1"     "BECN2"     "BIRC6"    
 [9] "BLOC1S1"   "BLOC1S2"   "BORCS5"    "BORCS6"    "BORCS7"    "BORCS8"    "CACNA1A"   "CALCOCO2" 
 ...
 
```

&nbsp;

## 3 Sample Data Analysis

```R
 
mySys <- getSysInteractions(sysName = "SLIGR", criterion = "stringent")

hypothesize(mySys)
```
![ppi_ggi](https://github.com/NElnour/BCB420.2019.ESA/blob/master/inst/extdata/exosc6.png?raw=true)
=======
## 3 Functions


```R
mySys2 <- getSysInteractions(sysName = "SLIGR", criterion = "relaxed")
hypothesize(mySys2, mySys)
```
![all_ggi](https://github.com/NElnour/BCB420.2019.ESA/blob/master/inst/extdata/network.png?raw=true)
&nbsp;

To look at several systems of interest at a time,
```R
systems <- c("PHALY", "SLIGR", "NLRIN", "HVGCR")
mySys <- getSysInteractions(systems)
mySys2 <- getSysInteractions(systems, criterion = "relaxed")

hypothesize(mySys)
hypothesize(mySys2, mySys)
```

![.](https://github.com/NElnour/BCB420.2019.ESA/blob/master/inst/extdata/multiSystemPPIGGI.png?raw=true)
![..](https://github.com/NElnour/BCB420.2019.ESA/blob/master/inst/extdata/multiSystemAllGGI.png?raw=true)
#### 3.1 Fetching data

See:
```R
?fetchData 
fetchData()
```
&nbsp;

#### 3.2 Systems database utilities


'./R/SyDButils.R' contains a number of exported functions to support work with a systems database:

* `SyDBgetRootSysIDs()` - returns a named vector with the IDs of all root systems in a systems database;
* `SyDBgetSysSymbols()` - returns a list of the same length as the input vector with HGNC symbols of each system in the input vector;
* `SyDBgetIDforKey()` - returns a vector of IDs for key(s) in a column of a table;
* `SyDBgetValforID()` - returns a vector of values in a column of a table where the ID matches the input;
* `SyDBinvalid()` - checks whether its database argument is a valid, current, system database;
* `SyDBTree()` - returns a tree representation of the hierarchical structure of a system or systems in the database.

Examples:
```R
mySDB <- fetchData("SysDB")

SyDBgetRootSysIDs(mySDB)
#                                       PHALY                                       SLIGR 
# "give.jams-1d8-648b-1e12-6a9f-65421424affe" "cast.rear-5f1-56ef-1cd2-1ae3-54a5556d59ff" 
#                                       NLRIN                                       HVGCR 
# "scar.blur-9bc-29cf-31f2-1981-4d92edf4d0e6" "help.mink-f96-e98b-9e12-ab41-8217a3ecb0cd" 


names(SyDBgetRootSysIDs(mySDB))
# [1] "PHALY" "SLIGR" "NLRIN" "HVGCR"


SyDBgetIDforKey("HOPS complex", "code", "component", mySDB)
# [1] "flip.face-782-d729-9622-299c-073fb8f2ec94"


SyDBgetSysSymbols(mySDB, "HOPS complex")
# $`HOPS complex`
# [1] "VPS11"  "VPS16"  "VPS18"  "VPS33A" "VPS39"  "VPS41" 


cat(SyDBTree("NLRIN", mySDB, MAX = 3), sep = "\n")
# 
#   --NLRIN
#     |___NLRP3 activation signal
#         |___Common activating signal
#         |___NLRP3
#         |___NLRP3 activation
#     |___NLRIN regulation
#         |___NLRIN positive regulation
#         |___NLRIN negative regulation
#     |___NLRP3 priming signal
#         |___IL-beta priming signal
#         |___TNF priming signal
#         |___TLR priming signal
#     |___NOT IN NLRIN
#         |___NLRP7
#         |___NLRC4
#         |___AIM2
#         |___NLRP14
#         |___NLRP6
#         |___NLRP12
#         |___NLRP2
#         |___NLRP9b
#         |___NLRP1b


```

&nbsp;


## 4 Notes

... in progress
&nbsp;

>>>>>>> upstream/master
## 4 References and Further Reading

... in progress
&nbsp;

## 5 Acknowledgements

... in progress
&nbsp;

<!-- [END] -->
