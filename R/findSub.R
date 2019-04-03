# findSub.R


#----prepare----
result_table <- matrix("",0,3)
colnames(result_table) <- c("System_Name", "New_Subsystem_Name", "Components")
map_complex <- readRDS("~/GitHub/BCB420.2019.ESA/data/map_complex.rds")
map_enrichment <- readRDS("~/GitHub/BCB420.2019.ESA/data/map_enrichment.rds")
map_prob_complex <- readRDS("~/GitHub/BCB420.2019.ESA/data/map_prob_complex.rds")
#---------------

#' \code{findSub} Find new subsystem for the system we asked.
#'
#' @param systemName (character) The system name that have 6 characters. All characters
#'                                       are uppercase.e.g. PHALLY
#' @param threshold (integer) The thershold for selection possible subsystem.
#' @param db (character) The database.
#'
#' @return (table) a table with possible subsystem name and its components.
#' @examples
#' findsub("PHALLY", 0.2)
#' @export
findSub <- function(systemName, threshold, db){
    result_table <- rbind(result_table, findSubInCom(systemName, threshold, db))
    allsubsystem <- fetchSubsystem(db, systemName)
    #loop all substem in a system
    if (length(allsubsystem) != 0){
    for (subsystem in allsubsystem){
      findSub(subsystem, threshold,db)
      }
    }
    return (result_table)
    }

#' \code{findSubInCom} Find subsystems in the provided components by using hu.map database.
#'
#' @param systemName (character) The system name that have 6 characters. All characters
#' @param db (character) The database.
#' @param threshold (integer) The thershold for selection possible subsystem
#' @examples
#' mySDB <- fetchData("SysDB")
#' findSubIncom("PHALLY", mySDB)
#' @return (character) A list of components in the system but not in the subsystem
#' @export
findSubInCom <- function(systemName, threshold, db) {
 # New subsystem must appear in compnents that don't belong to any other subsystem
 comNotinSub <- fectchComponentNotInSub(systemName, db)
 if (length(comNotinSub) ==0){
   return (NULL)
 }
 # Create possibility set for creating new subsystem
 possibleSubs <- c()
 finalPossibleSubs <- c()
 # List all enrichments of components that not in subsystem
 enrichmentCom <-c()
 for (components in comNotinSub){
   for (component in components){
   enrichment <- c()
   # Find complex it belongs to
   possibleCom <- which(component == map_complex)
   # find enrichments for each components in hu.map
   for (complex in possibleCom){
        enrichment <- union(enrichment,map_enrichment[which(toString(complex) == map_enrichment[,1]), 2])
   }
    enrichmentCom <- append(enrichmentCom,enrichment)
  }}
 # find repreated component and let them bind together and put them into the possible subsystem
 allEnrichment <- c()
 for (enrichments in enrichmentCom){
   list_enrichments <- strsplit(enrichments, "//")
   for (enrichment in list_enrichments){
   allEnrichment<- append(allEnrichment, enrichment)
 }}
 freqTable <- table(allEnrichment)
 filtedFreqTable <- freqTable[freqTable > 1]
 enrichmentNames <- names(sort(filtedFreqTable,decreasing=TRUE))

 for (name in enrichmentNames){
   possibleSub <- c()
   # find the position of component with the enrichment we asked
   for (i in 1:length(enrichmentCom)){
     list_enrichments <- strsplit(enrichmentCom[i], "//")
     for (enrichment in list_enrichments){
     if (name %in% enrichment){
        possibleSub <- paste(possibleSub,comNotinSub[i])
     }
     }
     }
   possibleSubs <- append(possibleSubs, possibleSub)
   }


 # find bi-complex possibility score from hu.map and calculate if the possibility score \
 # in a possibile subsystem is larger than theroshold]
 for (possibleSub in possibleSubs){
   list_possibleSub <- strsplit(possibleSub, " ")
   numCom <- length(list_possibleSub)
   for (i in 1:numCom){
     # find bi-complex possibility scores for each combination in possible subsystem
     for (n in (i+1):numCom){
       locations <- which(map_prob_complex[2,] == possibleSub[i], map_prob_complex[4,] == possibleSub[n])
       locations <- union(locations,which(map_prob_complex[4,] == possibleSub[i], map_prob_complex[2,] == possibleSub[n]))
       for (location in locations){
          if(strtoi(location[5], base = 10L) > threshold){
            break
          }
       }
     }

   }
   finalPossibleSubs <- append(finalPossibleSubs,possibleSub)
 }
 finalPossibleSubsNames <- c()
 # Find subsystem names for subsystem
 for (possibleSub in finalPossibleSubs){
   for (i in 1:length(possibleSubs)){
     if (possibleSub == possibleSubs[i]){
       finalPossibleSubsNames <- append(finalPossibleSubsNames,enrichmentNames[i])
     }
   }
 }
 # Return a table with system name, subsystem Names and possible components.
 if (length(finalPossibleSubsNames) != 0){
 result <- matrix("",0,3)
 colnames(result) <- c("System_Name", "New_Subsystem_Name", "Components")
 for (i in 1:length(finalPossibleSubs)){
 result <- rbind(result, c(systemName, finalPossibleSubsNames[i], finalPossibleSubs[i]))
 }
 return (result)
 } else {
   return (NULL)
 }
}

#' \code{fectchComponentNotInSub} Find components that don't appear in any subsystem
#'                                in the provided system .
#'
#' @param systemName (character) The system name that have 6 characters. All characters
#' @param db (character) The database.
#'
#' @examples
#' mySDB <- fetchData("SysDB")
#' fectchComponentNotInSub("PHALLY", mySDB)
#' @return (character) A list of components in the system but not in the subsystem
#' @export
fectchComponentNotInSub <- function(systemName, db) {
  allcom <- SyDBgetSysSymbols(db, systemName)
  # Get ride of old format
  for (allcom_ in allcom){
    allcom <- allcom_
  }
  allsubsystem <- fetchSubsystem(db, systemName)
  allComInsub <- c()
  for (subsystem in allsubsystem) {
      comSubs <- SyDBgetSysSymbols(db, subsystem)
      # Get ride of old format
      for (comSubs_ in comSubs){
        comSubs <- comSubs_
      }
      for (com in comSubs){
      allComInsub <- append(allComInsub,com)
    }}
  comNotInsub <- allcom[is.na(pmatch(allcom,allComInsub))]
  return (comNotInsub)
}

#' \code{fetchSubsystem} Find subsystem in the provided system.
#'
#' @param systemName (character) The system name that have 6 characters. All characters
#' @param db (character) The database.
#' @return (character) A list of subsystem in the system
#' @examples
#' mySDB <- fetchData("SysDB")
#' fetchSubsystem(mySDB, "PHALLY")
#'
#' @export
fetchSubsystem <- function(db, sys){
  tree <- SyDBTree(sys, db, MAX = 3)
  subsystem <- c()
  for (component in tree){
    if (startsWith(component, "    |___")){
      component <- substring(component,9,100000)
      subsystem <- union(subsystem, component)
  }}
  return (subsystem)
}

# [END]
