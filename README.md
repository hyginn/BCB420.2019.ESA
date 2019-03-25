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
###### [Nada Elnour](https://orcid.org/0000-0001-6165-1542),
###### University of Toronto, Canada
###### &lt; nada.elnour@mail.utoronto.ca &gt;

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
1. About this Package<br/>
2. Data<br/>
3. Notes<br/>
4. References<br/>
5. Acknowledgements<br/>
<!-- TOCabove -->

----

## 1. About this Package

This package is a joint development platform for student designed tools for Exploratory Systems Analysis, part of the University of Toronto course BCB420H1S (Computational Systems Biology) in the 2018/2019 academic year.

&nbsp;


----

## 2. Data

Suporting resources include curated systems data and other data resources:

#### 2.1 The HGNC symbol reference

Load the `HGNC` object into the data folder:

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


#### 2.3 STRING edges: Protein-protein Interaction Data

`STRINGedges` contains edges of the STRING database mapped to HGNC symbols.
Load `STRINGedges` into the data folder:

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

### 2.4 Expression profiles:

This is in progress. Here is a function stub that returns a random, repeatable, scaled expression profile for a HGNC gene symbol:

```R

exProf <- function(sym, hgnc = HGNC, ncol = 20) {
  # returns a set of numbers as a virtual expression profile, for 
  # development purposes only.
  set.seed(which(hgnc$sym == sym))
  p <- as.vector(scale(runif(ncol)))
  set.seed(NULL)
  return(p)
}

exProf("ARF5")
#   [1]  0.74593008  1.29119935  0.09176229  0.26580368  1.04807753
#   [6]  0.99590563 -1.58896376 -1.45692057  0.68341688 -0.09393910
#  [11]  0.72066664  0.60304597 -0.09448269  0.36606745 -1.58234271
#  [16]  0.30937106 -0.88034475 -1.33628536  1.23931446 -1.32728207

```

&nbsp;

#### 2.5 Systems:

This is in progress. Here is a function stub that returns a set of gene symbols for a system name:

```R

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

fetchComponents("PHALY")
#   [1] "AMBRA1"    "ATG14"     "ATP2A1"    "ATP2A2"    "ATP2A3"   
#   [6] "BECN1"     "BECN2"     "BIRC6"     "BLOC1S1"   "BLOC1S2"  
#  [11] "BORCS5"    "BORCS6"    "BORCS7"    "BORCS8"    "CACNA1A"  

fetchComponents("NONSUCH")
#  [1] ""

```

&nbsp;

#### 2.6 Sample Implementation --- Predicting Regulatory Networks from Undirected PPI
To begin working with a system, we need the excel sheet containing the systems components as parsed by [BCB420-2019-resources](https://github.com/hyginn/BCB420-2019-resources) scripts. For example, for the SLIGR system, the excel sheet is placed in the data folder:

```R
filename <- "../data/SLIGR.xlsx"
load("../data//HGNC.RData")

source('~/Documents/BCB420/BCB420.2019.ESA/R/predictDirectedNetworks.R', echo=FALSE) #enter credentials to generate your key

ensembl <- useMart(biomart = "ensembl")
human <- searchDatasets(mart = ensembl, pattern = "hsapiens")
myMart <- useMart("ensembl", human$dataset)

mySys <- getSysInteractions(filename, mart = myMart, criterion = "stringent")
head(mySys)

  gene1  gene2        interactionType
1   ATM   TP53 Phenotypic Enhancement
2  RELA CREBBP Phenotypic Enhancement
3  RELA CREBBP Phenotypic Suppression
4  RELA CREBBP Phenotypic Enhancement
5  RELA  EP300 Phenotypic Enhancement
6  RELA    FOS Phenotypic Enhancement
```
To get a hypothesized network, call hypothesize():
```R
hypothesize(mySys)
```
![PPI-GGI](https://github.com/NElnour/BCB420.2019.ESA/blob/master/inst/extdata/exosc6.png?raw=true)
&nbsp;

If instead we wanted both genetic and physical interactors in the SLIGR system,

```R
mySys2 <- getSysInteractions(filename, mart = myMart, criterion = "relaxed")

hypothesize(mySys2, mySys)
```
![relaxed_networks](https://github.com/NElnour/BCB420.2019.ESA/blob/master/inst/extdata/networks.png?raw=true)

&nbsp;

## 3. Notes

&nbsp;

## 4. References and Further Reading

&nbsp;

## 5. Acknowledgements

&nbsp;

<!-- [END] -->
