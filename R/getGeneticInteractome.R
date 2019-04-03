# getGeneticInteractome.R
#' @title getGeneticInteractome
#' @author Nada Elnour, \email{nada.elnour@@mail.utoronto.ca}
#'
#' \code{getGeneticInteractome} augments systems with BioGRID genetic
#' interaction annotations
#'
#' @param mySys (dataframe) A 2-column dataframe of PPI between system's
#' components.
#' @param criterion (string) Either "stringent" or "relaxed", specifying
#' whether only GGI between physical interactors should be selected.
#'
#' @return (dataframe) A 3-column dataframe of PPIs and either
#' the GGIs between them (iff \code{criterion} == "stringent") or
#' all GGIs that implicate them (iff \code{relaxed} == "stringent"). The first
#' two columns denote the interacting pair; the third is the type of genetic
#' interaction.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' biogrid <- fetchData("BioGRID")
#' url <- "http://www.informatics.jax.org/go/report.txt?goID=GO:0006099&results=42&startIndex=0&sort=term&dir="
#' myGenes <- fread(url)
#' myGenes <- myGenes$`MGI Gene/Marker ID`
#' myGenes <- toupper(myGenes)
#' PPI <- biogrid[biogrid$type == "physical", ]
#' sel <- PPI$A %in% myGenes & PPI$B %in% myGenes
#' mySys <- PPI[sel, ]
#' ppi_ggi <- getGeneticInteractome(mySys, "stringent")
#' network <- getGeneticInteractome(mySys, "relaxed")
#' hypothesize(network, ppi_ggi)
#' }
#' @export
getGeneticInteractome <- function(mySys, criterion){
  # load BioGRID fetched GGI
  file <- system.file("extdata", "humInt.txt",
                      package = "BCB420.2019.ESA",
                      mustWork = TRUE)
  humInt <- read.delim(file, stringsAsFactors = FALSE)
  humInt <- humInt[complete.cases(humInt), ]
  colnames(mySys)[1:2] <- c("a", "b")

  if (criterion == "stringent") { # stringent => GGI between physical
                                  # interactors of system components

    interctorsInSystemID <- (humInt$A %in% mySys$a & humInt$B %in% mySys$b)

  } else if (criterion == "relaxed") { # stringent => all GGI in which system's
                                       # physical interactors are involved
    interctorsInSystemID <- (humInt$A %in% mySys$a | humInt$B %in% mySys$b)
  }

  ggi <- humInt[interctorsInSystemID,]
  ggi <- ggi[complete.cases(ggi),]

  return(ggi)
}

# [END]
