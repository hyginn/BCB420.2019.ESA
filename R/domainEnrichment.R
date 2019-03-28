# domainEnrichment.R

#' Domain enrichment in a system.
#'
#' \code{domainEnrichment} Find the enriched domains in selected system and interpret enrichment by plot.
#'
#' @param system The shortehand notation of selected system, five-letter code as discussed in class. Default is PHALY.
#' @param numEnrichment The number of enriched domain ploted. Default is 40.
#' @return A data frame that include domains in the system with its counts and frequency of exist in system and in random genes; p-value, odd ratio, adjust p-value and enrichment; interpro descriptions.
#'
#' @import stats
#' @import biomaRt
#' @import ggplot2
#' @import ggiraph
#'
#' @author \href{https://orcid.org/0000-0001-8720-5874}{Fan Shen} (aut)
#'
#' @examples
#' # Plot 50 enrichend domains in PHALY system orderd by rich factor(counts in system/counts in random)
#' # User can hover around the point on plot to get more info (i.e counts, p-value)
#' system <- "PHALY"
#' domainEnrichment(system)
#' @export
domainEnrichment <- function(system = "PHALY", numEnrichment = 40) {
  # 1. Fetch data from Dr. Steipe's BCB420.2019.ESA package
  # Genes and their domains
  genesIPR <- fetchData("genesIPR")

  # domains and their genes
  IPRgenes <- fetchData("IPRgenes")

  # Get the genes for the system
  myDB <- fetchData("SysDB")
  genes <- SyDBgetSysSymbols(myDB, system)
  #genes <- fetchComponents(system)


  # 2. Find the enrichment counts
  # 2.1 store all the domains that in the system and their counts
  domains <- c()
  for (gene in genes) {
    domains <- c(genesIPR[gene][[1]], domains)
  }

  tmp <- as.data.frame(table(domains))
  colnames(tmp) <- c("domains", "sysCounts")

  # Store the contingency matrix for each gene
  contingency <- list()
  for(i in seq_along(tmp$domains)) {
    randomCount <- length(IPRgenes[as.character(tmp$domains[i])][[1]])
    tmp$randomCounts[i] <- randomCount
    contingency[[i]] <- c(tmp$sysCounts[i],
      randomCount,
      sum(genes %in% names(genesIPR))-tmp$sysCounts[i],
      length(genesIPR)-randomCount)
    names(contingency)[i] <- as.character(tmp$domains[i])
  }

  # 2.2 find the enrichment for each domain in the system
  tmp$sysFreq <- (tmp[,2]) / (sum(genes %in% names(genesIPR)))

  # 2.3 find the enrichment for each domain in all genes
  tmp$randomFreq <- tmp$randomCounts/length(genesIPR)


  # 3. Calculate the p_value for enrichment
  fisherTests <- data.frame(domain = character(),
                            p_value = numeric(),
                            odd_ratio = numeric(),
                            stringsAsFactors = FALSE)

  for (i in seq_along(contingency)) {
    fisherTest <- fisher.test(matrix(contingency[[i]],nrow=2,ncol=2),alternative="greater")
    add <- data.frame(domain = names(contingency)[i],
                      p_value = fisherTest$p.value,
                      odd_ratio = fisherTest$estimate,
                      row.names = i,
                      stringsAsFactors = FALSE)
    fisherTests[i,] <- add
  }

  fisherTests$p_adj <- stats::p.adjust(fisherTests$p_value)

  # Put all analyzed stats together
  tmp$enrichment <- (tmp$sysFreq)/(tmp$randomFreq)
  tmp <- merge(tmp, fisherTests, by.x = "domains", by.y = "domain")

  # 4. Find the description for enriched domains
  myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

  filters <- biomaRt::listFilters(myMart)
  attributes <- biomaRt::listAttributes(myMart)

  description <- biomaRt::getBM(filters = "interpro",
    attributes = c("interpro",
      "interpro_short_description",
      "interpro_description"),
    values = tmp$domains,
    mart = myMart)

  # Merge the two data frame, such that those with no description are excluded since they are family
  tmp <- merge(tmp, description, by.x = "domains", by.y = "interpro")

  # Reorder the tmp by enrichment
  tmp <- tmp[order(tmp$enrichment,decreasing = TRUE),]

  # Multiple testing?

  output <- tmp[1:numEnrichment,]
  output$tooltip <- c(paste0("Name = ", output$interpro_short_description,
                             "\n Counts in sys = ", output$sysCounts,
                             "\n Counts in random = ", output$randomCounts,
                             "\n p-value = ", output$p_value,
                             "\n odd ratio = ", output$odd_ratio,
                             "\n Enrichment = ", output$enrichment))
  # 5. Plot
  p1 <- ggplot(output, aes(x = output$sysCounts/output$randomCounts, y = reorder(output$interpro_short_description, output$sysCounts/output$randomCounts)))
  my_gg <- p1 +
    geom_point_interactive(tooltip = output$tooltip, size=2.5) +
    geom_point(aes(color=output$p_value)) +
    xlab("Rich Factor") +
    ylab("Interpro Description") +
    ggtitle(paste0("Enriched Domains in ", system))

  print(girafe(code = print(my_gg)))

  return(tmp)

  # TODO: Wrap up
  # describe the example a little bit!
  # use devtools::document() to get the namespace and md file
  # Add pkg for description
  # devtools::use_package("biomaRt") to add pkg
}

# [END]
