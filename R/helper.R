# helper.R
#'
#' \code{EnrichmentHelper} It's a helper function for pathway enrichment analysis 
#' function load all the required data for pathway enrichment analysis
#' Return a list that contain reactome information. Modified data annotation from 
#' reference: https://github.com/yoonsikp/BCB420.2019.REACTOME/
#'
#' @return A list of annotated information for r eactome
#' @author {Yufei Yang} (aut)
#' @examples
#' tmp<-EnrichmentHelper()
#'
EnrichmentHelper<-function(){
  HGNC <- fetchData("HGNCreference") 
  load(file = file.path("./data", "ensg2sym.RData"))
  #read reactome data
  tmp <- readr::read_delim(file.path("./data", "Ensembl2Reactome.txt"),
                           delim = "\t",
                           skip = 0,
                           col_names = c("ENSEMBL",
                                         "REACTOME_ID", 
                                         "HYPERLINK", 
                                         "Description", 
                                         "Unknown",
                                         "Species"))  
  tmp <- tmp[tmp$Species == "Homo sapiens",]
  tmp$HGNC <- ensg2sym[tmp$ENSEMBL]
  #remove some not used column
  drops <- c("REACTOME_ID","HGNC")
  tmp<-tmp[ , (names(tmp) %in% drops)]
  return(tmp)
}

