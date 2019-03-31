# getKey.R
#
# ==============================================================================
# SIDE EFFECTS:
# This script imports biogridr which uses the deprecated function xml2::xml_find_one().
# ==============================================================================
#
#' Generates a unique key for BioGrid API access
#'
#' @param my.name A string of the user's full name separated by space(s)
#' @param my.email A string matching a valid user email
#' @param my.project A string specifying a short project name (no spaces)
#' @return The unique key for user and project for data retrieval from BioGrid
#' @import biogridr
#' @examples
#' myKey <- getKey("Nada Elnour", "nada.elnour@@mail.utoronto.ca", "SILGRESA")
#'
#' @export
getKey <- function(my.name, my.email, my.project) {
  my.name <- unlist(strsplit(my.name, " "))

  options(warn = -1)
  myKey <-
    bg_get_key(my.name[1], my.name[length(my.name)], my.email, my.project)
  options(warn = 0)
  return(myKey)
}

# [END]
