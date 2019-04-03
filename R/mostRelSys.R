#mostRelSys.R
#
#purpose: find the most related system with the given gene name
#'Version:2.0.0
#'Data:2019.4
#'Author:Liwen Zhuang
#
#
#===PARAMETERS ====
#db: database which stores the given system data
#geneName: HGNC symbol of the interested gene
#
#
#===FUNCTIONS =====
#
#' RelSys
#'
#' Output dataframe with system names, number of conponents, number of interactions and score.
#'
#' @param db dataframe,database which stores the system data.
#' @param geneName string,the interested HGNC symbol
#' @return dataframe with information for each system: system name,number of components, number of interactions and score.
#' @export
#' @example
#'
#' STRINGedges <- fetchData("STRINGedges0.9")
#' myDB <- fetchData("SysDB")
#' RelSys(myDB,"CSK")
#'
RelSys <- function(db,geneName){
  numSys <- length(SyDBgetRootSysIDs(db))
  scoreindex <- 1
  SystemNameVec <- names(SyDBgetRootSysIDs(db))
  ComponentVec <- c()
  InteractionVec <- c()
  ScoreVec <- c()
  for (sysName in SystemNameVec){
    totalscore <- 0 #total score of edges between input and element in this system
    totaledges <- 0 #total number of edges existed
    SymbolList <- unlist(SyDBgetSysSymbols(db, sysName))#vectorize the list
    for (SysSymbols in SymbolList) {
      if (fetchScore(SysSymbols,geneName) != -1){#edge exist in normalized STRING database
        totalscore <- totalscore + fetchScore(SysSymbols,geneName)
        totaledges <- totaledges + 1
      }
    }
    finalScore <- totalscore/length(SymbolList)
    ComponentVec[scoreindex] <- length(SymbolList)
    InteractionVec[scoreindex] <- totaledges
    ScoreVec[scoreindex] <- finalScore
    scoreindex <- scoreindex + 1
  }
  res <- data.frame("SystemName"=SystemNameVec,
                  "Components"=ComponentVec,
                  "Interactions"=InteractionVec,
                  "Score"=ScoreVec)
  #order score from hight to low
  res <- res[order(res$Score,decreasing = TRUE),]
  return(res)
}


#' mostRelSys
#'
#' Output the name of system with highest score
#'
#' @param db dataframe,database which stores the system data.
#' @param geneName string,the interested HGNC symbol
#' @return String, name of system with highest score
#' @export
#' @example
#' STRINGedges <- fetchData("STRINGedges0.9")
#' myDB <- fetchData("SysDB")
#' mostRelSys(myDB,"CSK")
#'
mostRelSys <- function(db,geneName){
  df <- RelSys(db,geneName)
  return(toString(df[1,1]))
}


#'fetchScore
#'helper function for mostRelSys: fectch string score between g1 and g2
#'
#'@param g1 string,gene name of the first element
#'@param g2 string,gene name of the second element
#'@return double, confidence score between g1 and g2
#'
fetchScore <- function(g1,g2){
  STRINGedges <- fetchData("STRINGedges0.9")
  if (identical(STRINGedges$score[STRINGedges$b == g1 & STRINGedges$a == g2],numeric(0)) == FALSE){
    return(STRINGedges$score[STRINGedges$b==g1 & STRINGedges$a==g2])
  }
  else if (identical(STRINGedges$score[STRINGedges$b==g2 & STRINGedges$a==g1],numeric(0)) == FALSE){
    return(STRINGedges$score[STRINGedges$b==g2 & STRINGedges$a==g1])
  }
  else{# the edge does not exist in the STRING databse
    return(-1)
  }
}

#[END]
