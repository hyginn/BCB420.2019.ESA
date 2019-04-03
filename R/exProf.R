# exProf.R

#' Gene Expression Profile.
#'
#' \code{exProf} Returns the log-normalized expression values for gene provided.
#'
#'
#' @param sym (character) Symbol of gene being investigated
#' @param myQNXP (matrix) Matrix of expressions. Columns represent experiments and
#' rows individual genes
#'
#' @return (numeric) Vector of log-normalize expressions by condition
#'
#' @author \href{https://orcid.org/0000-0002-2821-6204}{Matthew McNei;} (aut)
#'
#' @export
#'
exProf <- function(sym, myQNXP = myQNXP) {
  # Returns a normalized gene expression profile

  # Select row of specific gene
  p <- myQNXP[sym,]
  # Names correspond to experiments
  experiments <- names(p)
  # Find the various cell line names, for comparison
  # We need to find each control for normalization
  # We would assume the control to be the nearest column which
  # contains ctrl or starts a new cell line
  cLines <- substr(experiments, 1, regexpr("\\.", experiments))
  # Find indexes where cell line changes
  controls <- which(cLines[-1] != cLines[-length(cLines)])
  # Add those indexes with a ctrl experiment
  controls <- c(1, controls + 1, which(regexpr("ctrl", experiments) > 0))
  controls <- sort(unique(controls))

  # cut up the columns by the nearest control experiment
  groups <- as.numeric(cut(1:length(p), breaks = controls, include.lowest = T, right = F))
  groups[is.na(groups)] <- length(controls)
  groups[controls[length(controls)]] <- length(controls)
  # Log normalize by the correct control
  normalized <- log(p/p[controls[groups]])
  # Remove the controls from the final conditions
  normalized <- normalized[-controls]
  return(normalized)
}

# [END]
