# installgenome.R
#

#' \code{installHumanGenomeAnnotation} install Human genome Annotation when necessary.
#'

#' @return (NULL)
#' #' @author {Yuhan Zhang} (aut)
#' @examples
#' \dontrun{
#' installHumanGenomeAnnotation()
#' }
#' @export
installHumanGenomeAnnotation <- function() {

    BiocManager::install("org.Hs.eg.db")

  return(invisible(NULL))
}

# [END]
