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
#' @param db (list)  a systems database
#'
#' @return a named vector of QQIDs; system codes are names
#' @examples
#' \dontrun{
#'SyDBgetRootSysIDs(mySDB)              # returns named IDs
#'names(SyDBgetRootSysIDs(mySDB))[2]    # returns the second system code
#' }
#' @export
SyDBgetRootSysIDs <- function(db) {
  # Return IDs of all systems that are not subsystems. IDs are a named
  # vector. Note: if a top-level system were to contain itself, it
  # would not be returned.
  if (as.logical(nchar(msg <- SyDBinvalid(db)))) { stop(msg) }

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
#' @param sys       (character)  system codes or system IDs
#' @param treatNOT  (character)  not implemented
#'
#' @return (list) a list of the same length as "sys", containing character
#'                vectors of HGNC symbols, or zero-length vectors if no
#'                symbols were found
#' @examples
#' \dontrun{
#' # all HGNC symbols of the PHALY system
#' SyDBgetSysSymbols(mySDB, "PHALY")
#'
#' # all genes in the database:
#' unlist(SyDBgetSysSymbols(mySDB, SyDBgetRootSysIDs(mySDB)), use.names = FALSE)
#' }
#' @export
SyDBgetSysSymbols <- function(db,
                              sys,
                              treatNOT = "exclude") {
  # Return all HGNC symbols contained as subsystems or components in system
  # with code "sysCode", recursively.
  #
  # ToDo:
  #    treatNOT: "exclude" - exclude "NOT IN ..." from output
  #    treatNOT: "include" - include "NOT IN ..." in output
  #    treatNOT: "split"   - include "NOT IN ..." in separate list elements
  #    treatNOT: "only"    - only list "NOT IN ..." component symbols
  #

  if (as.logical(nchar(msg <- SyDBinvalid(db)))) { stop(msg) }
  stopifnot(is.character(sys))

  getOneSysSymbols <- function(db, s) {

    QQIDpatt<-"^[a-z]{4}\\.[a-z]{4}-[0-9a-f]{3}-([0-9a-f]{4}-){3}[0-9a-f]{12}$"

    # process s as "code" or QQID:
    if (grepl(QQIDpatt, s)) {   # ... it's a QQID
      myQ <- s
    } else {                    # ... it should be a system code
      myQ <- SyDBgetIDforKey(s, "code", "system", db)
    }

    # fetch all components from the database
    sc <- db[["systemComponent"]][ , c("systemID", "componentID")]

    MAXSTEPS <- nrow(sc) + 1  # Counter to detect error condition. In case
    # this error is ever triggered, it points to an
    # error in the programming logic that could have
    # lead to an infinite loop - investigate.

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

    # fetch the gene's symbols
    sel <- db$gene$ID %in% myGenes & (! is.na(db$gene$symbol))
    mySym <- sort(unique(db$gene$symbol[sel]))

    return(mySym)
  }

  symList <- list()
  for (s in sys) {
    symList[[s]] <- getOneSysSymbols(db, s)
  }

  return(symList)
}



#' \code{SyDBgetIDforKey} fetch an ID from a systems database.
#'
#' \code{SyDBgetIDforKey} documentation in progress ...
#'
#' @param val       (character)  match value "val" ...
#' @param att       (character)  ... in attribute column "att" ...
#' @param tbl       (character)  ... of table "tbl" ...
#' @param db        (character)  in database "db".
#'
#' @return (character) a vector of HGNC symbols or a zero-length vector if
#'                     none were found
#' @examples
#' \dontrun{
#'SyDBgetIDforKey("PHALY","code","system", mySDB) # get the ID for the system
#'                                                # with code "PHALY"
#' }
#' @export
SyDBgetIDforKey <- function(val, att, tbl, db) {
  # Return the ID(s) that match value "val" of attribute "att" in table "tbl"
  # of system database "db", version "v".

  if (as.logical(nchar(msg <- SyDBinvalid(db)))) { stop(msg) }
  stopifnot(is.character(c(val, att, tbl)))

  sel <- which(db[[tbl]][ , att] == val)
  return(db[[tbl]][sel, "ID"])
}


#' \code{SyDBinvalid} validate a systems database.
#'
#' \code{SyDBinvalid} documentation in progress. Tests correct type, structure,
#'                    and version of the database ...
#'
#' @param db        (character)  s systems database "db".
#'
#' @return (character) an error message describing the problem if the database
#'                     is invalid, or an empty string if the database is valid.
#' @examples
#' \dontrun{
#' # use the following idiom in functions:
#' if (as.logical(nchar(msg <- SyDBinvalid(<database parameter>)))) {stop(msg)}
#' }
#' @export
SyDBinvalid <- function(db) {
  ob <- deparse(substitute(db))
  v <- "2.1.1"
  SyDBtables <- c("parameter", "type", "system", "systemComponent",
                  "component", "componentMolecule", "molecule",
                  "geneProduct", "gene", "note")
  if (is.null(db) || ! is.list(db)) {
    msg <- sprintf("\"%s\" is not a list.", ob)
  } else if ( ! setequal(names(db), SyDBtables)) {
    msg <- sprintf("\"%s\" does not have %s.",
                   ob, "the expected structure of a systems database")
  } else if ((x <- db$parameter[db$parameter$typeID ==
                                db$type[db$type$name == "DBversion", "ID"],
                                "value"])  != v) {
    msg <- sprintf("Systems db \"%s\" has version %s - we need %s.", ob, x, v)
  } else {
    msg <- ""
  }
  return(msg)
}


# [END]
