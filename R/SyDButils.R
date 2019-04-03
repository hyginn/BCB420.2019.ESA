# SyDButils.R
#
# Utility functions for work with a systems database
#   - SyDBgetRootSysIDs()
#   - SyDBgetSysSymbols()
#   - SyDBgetIDforKey()
#   - SyDBgetValForID
#   - SyDBinvalid()
#   - SyDBTree()

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
#' @param val       (character)  IDs matching value(s) "val" ...
#' @param att       (character)  ... in any attribute column "att" ...
#' @param tbl       (character)  ... of the table "tbl" ...
#' @param db        (character)  in database "db".
#' @param all       (logical)  if TRUE, then NA will be added to the returned
#'                             vector, one for each element in "val" that was
#'                             not found in "att". If FALSE (default), only
#'                             matching IDs will be returned, i.e. the length
#'                             of the input is not preserved.
#'
#' @return (character) a vector of QQIDs, or NA where
#'                     none were found
#' @examples
#' \dontrun{
#'# get the IDs for the systems with code "PHALY" and "nosuch". Element 2 is NA.
#'SyDBgetIDforKey(c("PHALY", "nosuch"), "code","system", mySDB, all = TRUE)
#' }
#' @export
SyDBgetIDforKey <- function(val, att, tbl, db, all = FALSE) {

  if (as.logical(nchar(msg <- SyDBinvalid(db)))) { stop(msg) }
  stopifnot(is.character(c(val, att, tbl)) && is.logical(all))
  stopifnot(length(tbl) == 1)

  v <- character()
  for (i in seq_along(val)) { # for every requested value
    sel <- apply(db[[tbl]][ , att, drop = FALSE],  # TRUE if any value in
                 1,                                  # att is equal to val[i]
                 FUN = function(x){any(x == val[i])})
    if (all == TRUE) {       # if input length has to be preserved ...
      if (sum(sel) > 1) {    # ... there can't be multiple matches
        stop(sprintf("Value %s is not unique in requested columns",
                     as.character(val[i])))
      } else {               # if we have zero or one match, get NA or the ID
        v[i] <- ifelse(sum(sel) == 0, NA, db[[tbl]][sel, "ID"])
      }
    } else { # ... input length does not need to be preserved
      v <- c(v, unique(db[[tbl]][sel, "ID"]))  # ... append all unique matched IDs
    }
  }
  return(v)
}


#' \code{SyDBgetValForID} fetch an an attribute value for an ID from a systems database.
#'
#' \code{SyDBgetValFforID} documentation in progress ...
#'
#' @param id        (character)  return values matching IDs in "id" ...
#' @param att       (character)  ... from column "att" ...
#' @param tbl       (character)  ... of the table "tbl" ...
#' @param db        (character)  in database "db".
#' @param all       (logical)  if TRUE, then NA will be added to the returned
#'                             vector, one for each element in "id" that was
#'                             not found. If FALSE (default), only values for
#'                             matching IDs will be returned, i.e. the length
#'                             of the input is not preserved.
#'
#' @return (character) a vector of values, or NA where
#'                     none were found
#' @examples
#' \dontrun{
#'# get the IDs for the systems with code "notAnID" and "PHALY". Element 1 is NA.
#'myIDs <- c("notAnID", "give.jams-1d8-648b-1e12-6a9f-65421424affe")
#'SyDBgetValForID(myIDs, "code", "system", mySDB, all = TRUE)
#' }
#' @export
SyDBgetValForID <- function(id, att, tbl, db, all = FALSE) {

  if (as.logical(nchar(msg <- SyDBinvalid(db)))) { stop(msg) }
  stopifnot(is.character(c(id, att, tbl)) && is.logical(all))
  stopifnot(length(tbl) == 1)  # exactly on column must exist

  v <- character()
  for (i in seq_along(id)) { # for every requested id
    sel <- db[[tbl]][ , "ID"] == id[i]
    if (all == TRUE) {       # if input length has to be preserved ...
      if (sum(sel) > 1) {    # ... there can't be multiple matches
        stop(sprintf("ID %s is not unique in requested table", id[i]))
      } else {               # if we have zero or one match, get NA or the ID
        v[i] <- ifelse(sum(sel) == 0, NA, db[[tbl]][sel, att])
      }
    } else { # ... input length does not need to be preserved
      v <- unique(c(v, db[[tbl]][sel, att]))  # ... append matched values
    }
  }
  return(v)
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


#' \code{SyDBTree} print a (sub)system as a tree
#'
#' \code{SyDBTree} Produce a tree that shows a (sub)system's hierarchical
#' composition. (Sub)systems that are not found in "db" are not shown. If
#' none of the requested input (sub)systems can be found, the function
#' returns NULL. Subsystems whose code begins with "NOT IN " are sorted to
#' the end of the tree.
#'
#' @param sys       (character)  vector of codes or QQIDs of systems
#' @param db        (character)  system database
#' @param MAX       (numeric)    the maximum number of levels to return.
#'                                 Defaults to Inf.
#'
#' @return (character) a vector that can be cat()ed to draw the tree
#'
#' @examples
#' \dontrun{
#'# use cat(..., sep = "\n") to print a tree for a subsystem
#' cat(SyDBTree(c("NF-kappa-B pathways"), mySDB), sep = "\n")
#' }
#' @export
SyDBTree <- function(sys, db, MAX = Inf) {

  if (as.logical(nchar(msg <- SyDBinvalid(db)))) { stop(msg) }
  stopifnot(is.character(c(sys)))
  stopifnot(is.numeric(c(MAX)))

  # function to create typographic output
  makeLines <- function(x) {
    # replace components of x with graphical elements for printing
    # a "tree", for all but the last component.
    l <- length(x)
    if (l == 1) {
      x <- c("\n  --", x)
    } else {
      x <- c(rep("    ", l - 1), "|___", x[l])
    }
    return(paste0(x, collapse = ""))
  }

  # process sys as "code" or QQID: "codes" can come from the "system" table or
  # from the "component" table.
  QQIDpatt<-"^[a-z]{4}\\.[a-z]{4}-[0-9a-f]{3}-([0-9a-f]{4}-){3}[0-9a-f]{12}$"
  myQ <- character()
  sel <- ( ! grepl(QQIDpatt, sys))
  myQ <- c(myQ, SyDBgetIDforKey(sys[sel], "code", "system",    db))
  myQ <- c(myQ, SyDBgetIDforKey(sys[sel], "code", "component", db))
  myQ <- unique(myQ[ ! is.na(myQ)])

  # working subset of "db"
  mySC <- db$systemComponent[, c("systemID", "componentID")]

  # Depth-first search:
  # "H" is a list of vectors that hold a leaf of a tree and its ancestry.
  # Start with H <- sys. Then move along iH list-pointer. At every H[[iH]],
  # append children. If there is more than one child, apend H[[iH]] to end of
  # list. Terminate when pointer runs past end of list, i.e. iH > length(iH).

  H <- as.list(myQ)
  iH <- 1

  MAXIH <- nrow(mySC) + length(myQ)  # Counter to detect error condition.
  # In case iH ever exceeds MAXIH, this points to an error in the programming
  # logic that could have lead to an infinite loop, or a circular systems
  # definition - investigate.

  while (iH <= length(H) && iH <= MAXIH) {
    # fetch parent node from list: it is the last element at the current
    # list pointer
    thisP <- H[[iH]][length( H[[iH]] )]

    # apend all its child-nodes to copies of the parent node; append to list
    sel <- mySC$systemID == thisP
    for (C in mySC$componentID[sel]) {
      H[[length(H) + 1]] <- c(H[[iH]], C)
    }
    iH <- iH + 1                              # set pointer to next list item
  }
  stopifnot(iH < MAXIH) # "i.e. Infinite loop? Circular system definition?"

  # restrict depth
  sel <- as.logical(sapply(H,
                           FUN <- function(x) {
                             ifelse(length(x) > MAX, FALSE, TRUE)
                           }))
  H <- H[sel]

  # To sort "NOT IN..." nodes to the end, first we prefix "z-" as a sort-to-end
  # token to all QQIDs that have a "NOT IN" code in the system table ...
  notInIDs <- db$system$ID[grepl("^NOT IN ", db$system$code)]
  H <- lapply(H, FUN = function(x) {sel <- x %in% notInIDs
  x[sel] <- paste0("z-", x[sel])
  return(x)})

  # sort
  if (length(H) > 1) {
    vO <- order(unlist(lapply(H, FUN = function(x) paste(x, collapse = ""))))
    H <- H[vO]
  }

  # ToDo: move "NOT IN ..." nodes to end

  # ToDo: switch the type of output here
  #
  # Output mode "code":
  myCol <- "code"

  # replace QQIDs with desired values
  for (i in seq_along(H)) {
    for (j in seq_along(H[[i]])) {
      id <- gsub("^z-", "", H[[i]][j]) # remove sort-to-end token from key
      s <- SyDBgetValForID(id = id,    # replace ID in components
                           att = myCol,
                           tbl = "component",
                           db = db)
      if (length(s) == 0) {
        s <- SyDBgetValForID(id = id,  # replace ID in root systems
                             att = myCol,
                             tbl = "system",
                             db = db)
      }
      if (length(s) == 0) {            # replace ID for unknown IDs
        s <- "???"
      }
      H[[i]][j] <- s
    }
  }

  # replace components with graphical elements
  H <- unlist(lapply(H, makeLines))

  return(H)
}


# [END]
