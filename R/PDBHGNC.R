#' PDBHGNC.R
#'
#' The PDB-HGNC ID mapping data
#'
#' This is a dataframe holding the PDB-HGNC annotated database from the
#' BCB420.2019.PDB package from a previous BCB420 project by Heewon Lee.
#'
#' @format a data frame contianing HGNC and PDB transcript mappings
#' \describe{
#'   \item{HGNC}{HGNC entries mapped to transcript IDs}
#'   \item{Transcripts}{HGNC transcripts mapped to PDB chain IDs}
#'   \item{pdbChains}{PDB chain IDs mapped to HGNC transcripts and PDB IDs}
#'   \item{PDB}{PDB IDs mapped to PDB chain IDs}
#' }
#' @source \url{https://github.com/judyheewonlee/BCB420.2019.PDB}
#' @docType data
#' @name PDBHGNC
NULL

# [END]

