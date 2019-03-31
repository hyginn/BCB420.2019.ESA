# fetchPDBHGNCdatabase.R
#
# Purpose: Helper function in Exploratory System Analysis tool
# project for BCB420 2019
# Version: 1.0.0
# Version history: 1.0.0
# Date: March 28 2019
# Author: Heewon Lee (heewon.lee@mail.utoronto.ca)
# License: MIT
#
# Output: A dataframe contianing the PDB-HGNC annotated database
# linking PDB ID to HGNC ID with annotations from the
# BCB420.2019.PDB package. This is used in the SeqComparisonTable() function
# in the BCB420.2019.ESA package.
#
# Dependencies: BCB420.2019.PDB, BiocManager, msa, data.table, xml2, rtracklayer,
# biomaRt, devtools

# ====  FUNCTION  =============================================================
# fetchPDBHGNCdatabase.R
#'
#' \code{fetchPDBHGNCdatabase} Retrieve the database containing
#' the annotated database for HGNC and PDB relations created by author
#' Heewon Lee from previous BCB420 course project in the
#' BCB420.2019.PDB package on GitHub linked here:
#' https://github.com/judyheewonlee/BCB420.2019.PDB
#'
#' @return A dataframe contianing the PDB-HGNC annotated database
#' linking PDB ID to HGNC ID with annotations from the
#' BCB420.2019.PDB package. This is used in the SeqComparisonTable() function
#' in the BCB420.2019.ESA package.
#'
#' @author {Heewon Lee} (aut)
#' @examples
#'
#' # Call the fetchPDBHGNCdatabase function to generate a dataframe
#' # containing the annotated PDB-HGNC database from the
#' # BCB420.2019.PDB package
#' pdbHGNCdata <- fetchPDBHGNCdatabase()
#'
#' @export

fetchPDBHGNCdatabase <- function() {

  ### =========== Load all required packages and data ================= ###

  if (! requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }

  devtools::install_github("judyheewonlee/BCB420.2019.PDB", quiet = TRUE)

  if (! requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }

  if (! requireNamespace("xml2", quietly = TRUE)) {
    install.packages("xml2")
  }

  if (! requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  if (! requireNamespace("rtracklayer", quietly = TRUE)) {
    BiocManager::install("rtracklayer")
  }

  if (! requireNamespace("biomaRt", quietly = TRUE)) {
    BiocManager::install("biomaRt")
  }

  # Call PDBdataset() function from the BCB420.2019.PDB
  # package to assign pdbHGNC as the dataframe containing
  # the annotated database

  pdbHGNC <- BCB420.2019.PDB::PDBdataset()

  # Return the generated database
  return(pdbHGNC)

}



# [END]

