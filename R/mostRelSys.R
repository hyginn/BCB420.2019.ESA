#mostRelSys.R
#
#purpose: find the most related system with the given gene name
#Version:1.0.0
#Data:2019.3
#Author:Liwen Zhuang
#
#Input:System data, normalized STRING edge score, and gene name of interest
#Output:The most related system in the given system data
#Dependencies:None
#
#
#===PARAMETERS ====
#db: database which stores the given system data
#geneName: HGNC symbol of the interested gene
#
#
#===FUNCTIONS =====
#
#' mostRelSys
#'
#' main funcion, output the most related system of a given gene.
#'
#' @param db dataframe,database which stores the system data.
#' @param geneName string,the interested HGNC symbol
#' @return A name of the System which is most related to the given gene name.
#' @export
#' @example
#'
#' myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/","STRINGedges-2019-03-14.RData")
#' load(url(myURL))
#'
#' myDB <- fetchData("SysDB")
#' mostRelSys(myDB,"CSK")
#'
mostRelSys <- function(db,geneName){
  numSys<- length(names(SyDBgetRootSysIDs(db)))
  scoreindex<-1
  scorelist<-c()
  for (sysName in names(SyDBgetRootSysIDs(db))){
    totalscore<-0 #total score of edges between input and element in this system
    totaledges<-0 #total number of edges existed
    SymbolList<-unlist(SyDBgetSysSymbols(db, sysName))#vectorize the list
    for (SysSymbols in SymbolList) {
      if (fetchScore(SysSymbols,geneName) != -1){#edge exist in normalized STRING database
        totalscore<-totalscore+ fetchScore(SysSymbols,geneName)
        totaledges<-totaledges+1
      }
    }
    finalScore <- totalscore/length(SymbolList)
    scorelist[scoreindex]<-finalScore
    scoreindex<-scoreindex+1
    cat("total interaction found in system ",sysName,"is ",totaledges,"with final score of ",finalScore,".\n")
  }
  return(names(SyDBgetRootSysIDs(db))[which.max(scorelist)])
}


#'fetchScore
#'helper function for mostRelSys: fectch string score between g1 and g2
#'
#'@param g1 string,gene name of the first element
#'@param g2 string,gene name of the second element
#'@return double, confidence score between g1 and g2
#'
fetchScore <-function(g1,g2){
  if (identical(STRINGedges$score[STRINGedges$b == g1 & STRINGedges$a == g2],numeric(0))==FALSE){
    return(STRINGedges$score[STRINGedges$b==g1 & STRINGedges$a==g2])
  }
  else if (identical(STRINGedges$score[STRINGedges$b==g2 & STRINGedges$a==g1],numeric(0))==FALSE){
    return(STRINGedges$score[STRINGedges$b==g2 & STRINGedges$a==g1])
  }
  else{# the edge does not exist in the STRING databse
    return(-1)
  }
}

#====PROCEDURE===
#
#
#load data
myURL <- paste0("http://steipe.biochemistry.utoronto.ca/abc/assets/",
                "STRINGedges-2019-03-14.RData")
load(url(myURL))
myDB <- fetchData("SysDB")
#use the function
mostRelSys(myDB,"CFTR")
#> mostRelSys(myDB,"CFTR")
#total interaction found in system  PHALY is  4 with final score of  56.23077 .
#total interaction found in system  SLIGR is  0 with final score of  0 .
#total interaction found in system  NLRIN is  0 with final score of  0 .
#total interaction found in system  HVGCR is  6 with final score of  263.9048 .
#[1] "HVGCR"
