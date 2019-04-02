# domainEnrichment.R

#' Domain enrichment in a desired system.
#'
#' \code{domainEnrichment} Find the enriched domains in selected system and interpret enrichment by plot.
#'
#' For the PHYLA system, see details in the \href{../doc/domainEnrichment-Vignette.html}{PHALY domain enrichemnt vignette} (or load the vignette with \code{vignette("domainEnrichment-Vignette", package = "BCB420.2019.ESA")}).
#'
#' @section Details: This function is created as part of BCB420 2019 winter course. The domain databaase were fetched from InterPro, see more details in \code{\link{fetchData}}.
#'
#' @param system The shortehand notation of selected system, five-letter code as discussed in class.
#' @param alpha The alpha value for Benjamini-Hochberg control. Default is 0.05.
#' @return A data frame that include domains in the system with its counts and frequency in system and in random genes; p-value, odd ratio, adjust p-value and enrichment; interpro descriptions. Return NULL if system is not available by \code{\link{SyDBgetSysSymbols}}.
#'
#' @importFrom stats p.adjust fisher.test reorder
#' @import biomaRt
#' @import ggplot2
#' @import ggiraph
#'
#' @author \href{https://orcid.org/0000-0001-8720-5874}{Fan Shen} (aut)
#'
#' @seealso \code{\link{fetchData}} Fetch data from package
#' @seealso \code{\link{SyDBgetRootSysIDs}} Get the system ID in database
#' @seealso \code{\link{SyDBgetSysSymbols}} Get the genes in desired system
#'
#' @examples
#' # Plot enrichend domains in PHALY system orderd by rich factor
#' # Plot would separate to two parts if too many domains are enriched
#' # User can hover around the point on plot to get more info (i.e ID descripton)
#' system <- "PHALY"
#' domainEnrichment(system)
#'
#' @export
domainEnrichment <- function(system, alpha = 0.05) {
  # 1. Fetch data from Dr. Steipe's BCB420.2019.ESA package
  # Two lists from InterPro database
  genesIPR <- fetchData("genesIPR")
  IPRgenes <- fetchData("IPRgenes")

  # Genes for the system
  myDB <- fetchData("SysDB")
  if (! system %in% names(SyDBgetRootSysIDs(myDB))) {
    print("Input system should be in database of ESA package.")
    return(NULL)
  }

  genes <- SyDBgetSysSymbols(myDB, system)[[1]]


  # 2. Find the counts of genes for each domain
  # 2.1 Store all the domains that in the system and their counts
  domains <- c()
  for (gene in genes) {
    domains <- c(genesIPR[gene][[1]], domains)
  }

  tmp <- as.data.frame(table(domains))
  colnames(tmp) <- c("domains", "sysCounts")

  # 2.2 Store the counts from random (all genes)
  for(i in seq_along(tmp$domains)) {
    tmp$randomCounts[i] <- length(IPRgenes[as.character(tmp$domains[i])][[1]])
  }


  # 3. Data analysis
  # 3.1 Calculate the frequency for each domain in the system
  tmp$sysFreq <- (tmp[,2]) / (sum(genes %in% names(genesIPR)))

  # 3.2 Calculate the frequency for each domain in random (all genes)
  tmp$randomFreq <- tmp$randomCounts/length(genesIPR)

  # 3.3 Calculate the p_value and odd_ratio
  # Null hypothesis: The categories are independent.
  # Some useful page for understanding recommended by Dr. Steipe.
  # https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
  # https://www.pathwaycommons.org/guide/primers/statistics/multiple_testing/

  # Fisher Test Ref:
  # M. Berk (https://stats.stackexchange.com/users/30201/m-berk), Which
  # statistical test should be used to test for enrichment of gene lists?,
  # URL (version: 2013-10-11): https://stats.stackexchange.com/q/72556
  for (i in seq_along(tmp$domains)) {
    contingency <- c(tmp$sysCounts[i],
                     tmp$randomCounts[i],
                     sum(genes %in% names(genesIPR))-tmp$sysCounts[i],
                     length(genesIPR)-tmp$randomCounts[i])

    fisherTest <- fisher.test(matrix(contingency, nrow=2, ncol=2),
                                     alternative="greater")

    tmp$p_value[i] <- fisherTest$p.value
    tmp$odd_ratio[i] <- fisherTest$estimate
  }

  # 3.4 Multiple testing: Benjamini-Hochberg control for false discovery rate
  tmp$p_adj_BH <- p.adjust(tmp$p_value, "BH")

  # 3.5 Calculate the enrichment
  tmp$enrichment <- (tmp$sysFreq)/(tmp$randomFreq)


  # 4. Find the description for enriched domains
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  description <- biomaRt::getBM(filters = "interpro",
                                attributes = c("interpro",
                                               "interpro_short_description",
                                               "interpro_description"),
                                values = tmp$domains,
                                mart = myMart)

  # Merge the two data frame
  # note: domain IDs with no description are excluded since they are family
  tmp <- merge(tmp, description, by.x = "domains", by.y = "interpro")


  # 5. plot
  # 5.1 Get rid of false discovery.
  # Only tmp$p_adj_BH <= alpha. Small p-value means reject null hypothesis.
  output <- tmp[tmp$p_adj_BH <= alpha, ]

  # 5.2 Add output text for each domain
  output$tooltip <- c(paste0("ID = ", output$domains,
                             "\n Name = ", output$interpro_short_description,
                             "\n Counts in sys = ", output$sysCounts,
                             "\n Counts in random = ", output$randomCounts,
                             "\n p-value = ", output$p_value,
                             "\n odd ratio = ", output$odd_ratio,
                             "\n p-value BH = ", output$p_adj_BH,
                             "\n Enrichment = ", output$enrichment))
  output <- output[order(output$sysCounts/output$randomCounts, decreasing = TRUE), ]

  # 5.3 Plot the data
  if (nrow(output) > 40) {
    separate = ceiling(nrow(output)/2)
    for (i in seq(2)) {
      start <- (i-1)*separate+1
      output$label[start:nrow(output)] <- rep(paste0("Part ",i), nrow(output)-start+1)
    }
  } else {
    output$label <- rep("",nrow(output))
  }

  p <- ggplot(output,
               aes(x = output$sysCounts/output$randomCounts,
                   y = reorder(output$domains,
                               output$sysCounts/output$randomCounts)),scales = "free_y")
  my_gg <- p +
           geom_point_interactive(tooltip = output$tooltip, size=2.5) +
           geom_point(aes(color=output$p_value)) +
           labs(colour="P value")+
           xlab("Rich Factor") +
           ylab("Interpro ID") +
           ggtitle(paste0("Enriched Domains in ", system)) +
           theme(axis.text.y = element_text(angle = 20, hjust = 1)) +
           facet_wrap(~output$label, scales = "free")

  print(girafe(code = print(my_gg)))

  return(tmp)
}

# [END]
