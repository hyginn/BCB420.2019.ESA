# SyDButils.R
#
# Utility functions for work with a systems database
#   - SyDBgetRootSysIDs()
#   - SyDBgetSysSymbols()
#   - SyDBgetIDforKey()

#' \code{SyDBgetRootSysIDs} fetch root-system IDs from a systems database.
#'
#' \code{SyDBgetRootSysIDs} documentation in progress ..
#'
#' @param db (character)  name of the database
#' @param v  (character)  the database version. Defaults to 2.1.1
#'
#' @return a named vector of QQIDs; system codes are names
#' @examples
#' \dontrun{
#'SyDBgetRootSysIDs(mySDB)              # returns named IDs
#'names(SyDBgetRootSysIDs(mySDB))[2]    # returns the second system code
#' }
#' @export


SyDBgetRootSysIDs <- function(db, v = "2.1.1") {
  # Return IDs of all systems that are not subsystems. IDs are a named
  # vector. Note: if a top-level system were to contain itself, it
  # would not be returned.

  # get system IDs
  sIDs <- db[["system"]][ , "ID"]

  # which of those are not also a component?
  sIDs <- sIDs[ ! sIDs %in% db[["systemComponent"]][ , "componentID"]]

  sel <- db[["system"]][ , "ID"] %in% sIDs
  names(sIDs) <- db[["system"]][sel, "code"]
  return(sIDs)
}



#' \code{SyDBgetSysSymbols} fetch root-system IDs from a systems database.
#'
#' \code{SyDBgetSysSymbols} documentation in progress ...
#'
#' @param db        (character)  name of the database
#' @param sysCode   (character)  the requested system code
#' @param treatNOT  (character)  not implemented
#' @param v         (character)  the database version. Defaults to 2.1.1
#'
#' @return (character) a vector of HGNC symbols or a zero-length vector if
#'                     none were found
#' @examples
#' \dontrun{
#'SyDBgetSysSymbols("PHALY")            # all HGNC symbols in the PHALY system
#' }
#' @export

SyDBgetSysSymbols <- function(db,
                              sysCode,
                              treatNOT = "exclude",
                              v = "2.1.1") {
  # Return all HGNC symbols contained as subsystems or components in system
  # with code "sysCode", recursively.
  #
  # ToDo:
  #    treatNOT: "exclude" - exclude "NOT IN ..." from output
  #    treatNOT: "include" - include "NOT IN ..." in output
  #    treatNOT: "only"    - only list "NOT IN ..." component symbols
  #

  # fetch all components
  sysID <- SyDBgetIDforKey(sysCode, "code", "system", db)
  sc <- db[["systemComponent"]][ , c("systemID", "componentID")]

  MAXSTEPS <- nrow(sc) + 1  # Counter to detect error condition. In case
  # this error is ever triggered, it points to an
  # error in the programming logic that could have
  # lead to an infinite loop - investigate.

  myQ <- sysID  # queue that holds system IDs to process
  myComponents <- character()
  count <- 0
  while (length(myQ) > 0 && count <= MAXSTEPS) {
    count <- count + 1

    thisID <- myQ[1]                          # fetch No. 1 component from queue
    myComponents <- c(myComponents, thisID)   # add it to component list
    sel <- sc[ , "systemID"] == thisID        # add its sub-components to queue
    myQ <- c(myQ, sc[sel, "componentID"])
    sc <- sc[! sel, ]                         # remove processed components
    myQ <- myQ[-1]                            # remove No. 1 from queue
  }
  stopifnot(count <= MAXSTEPS)
  myComponents <- unique(myComponents)

  # which of these are atomic
  sel <- (db$component$ID %in% myComponents) &
    (db$component$componentType == "atomic")
  myComponents <- db$component$ID[sel]

  # collect all of their molecules
  sel <- db$componentMolecule$componentID %in% myComponents
  myMolecules <- db$componentMolecule$moleculeID[sel]

  # which of those are protein
  sel <- (db$molecule$ID %in% myMolecules) &
    (db$molecule$moleculeType == "protein")
  myProteins <- db$molecule$ID[sel]

  # collect all of their genes
  sel <- db$geneProduct$moleculeID %in% myProteins
  myGenes <- db$geneProduct$geneID[sel]

  # fetch their symbols
  sel <- db$gene$ID %in% myGenes & (! is.na(db$gene$symbol))
  mySym <- sort(unique(db$gene$symbol[sel]))

  return(mySym)
}



#' \code{SyDBgetIDforKey} fetch an ID from a systems database.
#'
#' \code{SyDBgetIDforKey} documentation in progress ...
#'
#' @param val       (character)  match value "val" ...
#' @param att       (character)  ... in attribute column "att" ...
#' @param tbl       (character)  ... of table "tbl" ...
#' @param db        (character)  in database "db".
#' @param v         (character)   database version. Defaults to 2.1.1
#'
#' @return (character) a vector of HGNC symbols or a zero-length vector if
#'                     none were found
#' @examples
#' \dontrun{
#'SyDBgetIDforKey("PHALY","code","system", mySDB) # get the ID for the system
#'                                                # with code "PHALY"
#' }
#' @export

SyDBgetIDforKey <- function(val, att, tbl, db, v = "2.1.1") {
  # Return the ID(s) that match value "val" of attribute "att" in table "tbl"
  # of system database "db", version "v".

  sel <- which(db[[tbl]][ , att] == val)
  return(db[[tbl]][sel, "ID"])
}



# [END]
