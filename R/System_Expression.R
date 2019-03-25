# System_Expression.R

#' System level expression values.
#'
#' \code{exProf} Returns the log-normalized expression values for gene provided.
#'
#' \code{System_Expression} Provides the mean expression level for each experiment
#' present in the dataset. As well, displays which conditions are up and down
#' regulted significantly.
#'
#'
#' @param sysname (character) 5-letter code of system being investigated
#' @param myQNXP (matrix) Matrix of expressions. Columns represent experiments and
#' rows individual genes
#'
#' @return (list) A list containing 3 elements. The first the mean log normalized
#' gene expression by condition. Second and third the conditions where the system
#' is significantly up and down regulated respectively.
#'
#' @author \href{https://orcid.org/0000-0002-2821-6204}{Matthew McNei;} (aut)
#'
#' @export

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

fetchComponents <- function(sys) {
  # returns a fixed set of symbols.
  # Function stub for development purposes only.
  if (sys == "PHALY") {
    s <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
           "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
           "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
           "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
           "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
           "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
           "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
           "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
           "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
           "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
           "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
           "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
           "VPS41", "VTI1B", "YKT6")
  } else {
    s <- ""
  }
  return(s)
}

System_Expression<- function(sysname, myQNXP){
  # Get a list of genes present in the given system
  components <- fetchComponents(sysname)

  # Normalized expression values for each gene in the system
  expressions <- sapply(components, exProf, myQNXP)

  # Find the means and variance of the expression changes caused by each condition
  meanExpression <- rowMeans(expressions, na.rm = T)
  varExpression <- apply(expressions, 1, var, na.rm = T)

  # calculate the one-sided p-values for each experiment
  p.values <- pnorm(meanExpression, sd = sqrt(varExpression))
  p.values[p.values > 0.5] <- 1 - p.values[p.values > 0.5]

  # Find the experiments where the system was significantly up
  # or down regulated at the 0.05 level
  UpConditions <- which(meanExpression > 0 & p.values < 0.05)
  DownConditions <- which(meanExpression < 0 & p.values < 0.05)

  experiments <- rownames(expressions)
  # Package results into list for export
  SysExpression <- list(`Mean Expression` = meanExpression,
                        `Upregulated Conditions` = experiments[UpConditions],
                        `Downregulated Conditions` = experiments[DownConditions])

  return(SysExpression)
}

# [END]
