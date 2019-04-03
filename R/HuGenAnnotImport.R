# HuGenAnnotImport.R

#' \code{HuGenAnnotImport} imports the org.Hs.eg.db, implemented to reduce import burden of package
#' This database will only be imported when the semanticSimilarity.R function is required
#'
#' @return (NULL) NULL
#' @examples
#' \dontrun{
#' HuGenAnnotImport()
#' }
#' @export

HuGenAnnotImport <- function() {
  if (! requireNamespace("org.Hs.eg.db")) {
    BiocManager::install("org.Hs.eg.db")
  }
  return(invisible(NULL))
}
