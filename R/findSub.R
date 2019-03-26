# makeSeq.R

#' \code{makeSeq} output a random DNA sequence.
#'
#' \code{makeSeq} outputs a string that contains a random sequence of characters
#' (A/C/G/T by default). By default the first three letters are "ATG" and the
#' last three letters are a stop-codon ("TAA", "TAG" or "TGA"). The number of
#' characters in the output is three times the \code{len} argument, i.e. it
#' contains \code{len} codons
#'
#' @param len (integer)       number of codons in the returned string
#' @param useInit  (logical)  make the first codon an initiation codon. Default
#'                            is \code{TRUE}.
#' @param useStop  (logical)  make the last codon a stop codon. Default is
#'                            \code{TRUE}.
#' @param alphabet (character vector)  the elements that are randomly
#'                                     concatenated. Default is A/C/G/T.
#' @param p (numeric vector)  Probabilites for each character to occur. Default
#'                            is equiprobable.
#' @param seed  (numeric)     the seed for the RNG. Default is \code{NULL} i.e.
#'                            the RNG behaves as if no seed had been set..

#' @return (character) a single string
#' @examples
#' makeSeq(7)
#' makeSeq(7, p = c(0.2, 0.4, 0.4, 0.2), seed = 112358)
#' gsub("T", "U", makeSeq(7)) # for RNA
#' @export
findSub <- function(systemName, threshold){
    curSystem <- systemName
    findSubInCom(curSystem, threshold)
    allsubsystem <- fetchSubsystem(systemName)
    #loop all substem in a system
    for (subsystem in allsubsystem){
      findSub(subsystem, threshold)
       }
}
findSubInCom <- function(systemName, threshold) {
 # New subsystem must appear in compnents that don't belong to any other subsystem
 comNotinSub <- fectchComponentNotInSub(systemName)
 # Create possibility set for creating new subsystem
 possibleSubs <- c()
 finalPossibleSubs <- c()
 # List all enrichments of components that not in subsystem
 enrichmentCom <-c()
  for (component in comNotinSub){
  # find enrichment in hu.map
    line <- grep(component,huMAPdataEnrichment)
   if (! is.null(line)){
     enrichment <- line[3]
     enrichmentCom <- append(enrichment, enrichmentCom)
   }
  }
 # find repreated component and let them bind together and put them into the possible subsystem
 freqTable <- table(enrichmentCom)
 filtedFreqTable <- freqTable[freqTabl > 1]
 enrichmentNames <- names(sort(filtedFreqTable,decreasing=TRUE))
 for (name in enrichmentNames){
  # find the position of component with the enrichment we asked
   locations <- match(enrichmentCom, name)
   possibleSub <- c()
   for (location in locations) {
     possibleSub <- append(possibleSub, comNotinSub)
   }
   possibleSubs <- append(possibleSubs, possibleSub)
 }
 # find bi-complex possibility score from hu.map and calculate if the possibility score \
 # in a possibile subsystem is larger than theroshold]
 for (possiblesub in possibleSubs){
   numCom <- length(possiblesub)
   for (i in 1:numCom){
     # find bi-complex possibility scores for each combination in possible subsystem
     for (n in (i+1):numCom){
       line <- which(huMapBiComplex[1,] == possibleSubs[i], huMapBiComplex[2,] == possibleSubs[n])
       if (line[3] <  threshold){
         break
       }
   }
   }
   finalPossibleSubs <- append(finalPossibleSubs,possiblesub)
 }
 # save finalPossibleSubs as new subsystems
}

fectchComponentNotInSub <- function(systemName) {
  allcom <- fetchComponents(systemName)
  allsubsystem <- fetchSubsystem(systemName)
  allComInsub <- c()
    for (subsystem in allsubsystem) {
      comSub <- fetchComponents(subsystem)
      allComInsub <- append(comSub, allComInsub)
    }
  comNotInsub <- allcom[is.na(pmatch(allcom,allComInsub))]
  return (comNotInsub)
}

# [END]
