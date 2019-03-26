# # Script to load STRING coexpression data and map it
# # v1 - 22/03/19
# # mapSTRING.R
# #
# #
# # A script that produces high-confidence, HGNC-symbol mapped datasets for both the
# # "9606.protein.links.detailed.v11.0.txt" and "9606.protein.actions.v11.0.txt"
# # data files available on the STRING consortium website
# #
# # Author: Boris Steipe (ORCID: 0000-0002-1134-6758)
# # Contributor: Gabriela Morgenshtern (ORCID:0000-0003-4762-8797)
# # License: (c) Author (2019) + MIT
# # Date: 2019-03-22
# #
# # ToDo:
# # Notes:
# #
# # ==============================================================================
#
# # SIDE EFFECTS: RData file creation
# # This script can be source()'d to create the RData files needed for default
# # package function.
# # All other code is executed when source()'d and this REQUIRES the user to
# # have the "9606.protein.links.detailed.v11.0.txt" and
# # "9606.protein.actions.v11.0.txt" in their /data directory.
#
# #### Most code comes from Dr. Steipe's BCB420.2019.STRING package: ####
# #### 1. Dependencies ####
# if (! requireNamespace("readr")) {
#   install.packages("readr")
# }
# if (! requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# if (! requireNamespace("biomaRt", quietly = TRUE)) {
#   BiocManager::install("biomaRt")
# }
#
# if (! requireNamespace("devtools")) {
#   install.packages("devtools")
# }
#
# # Question for B Steipe: why didn't this download the mapping tool script into /inst?
#
# if (! requireNamespace("devtools")) {
#   install.packages("devtools")
#   devtools::install_github("hyginn/BCB420.2019.STRING")
# }
#
# if (! requireNamespace("BiocCheck")) {
#   BiocManager::install("BiocCheck")
# }
#
#
#
# # Load raw STRING detailed dataset
# tmp <- readr::read_delim(file.path("./data", "9606.protein.links.detailed.v11.0.txt"),
#                          delim = " ",
#                          skip = 1,
#                          col_names = c("protein1", "protein2", "neighborhood", "fusion",
#                                        "cooccurence", "coexpression", "experimental", "database",
#                                        "textmining", "combined_score")) # 11,759,454 rows
# # Keeping col revlevant to our analysis (high coexpression edges)
# tmp <- tmp[,c("protein1", "protein2","coexpression","combined_score")]
# tmp <- tmp[(tmp$combined_score >= 800), ] #should combined_score be included?
#
# # Do all elements have the right tax id?
# all(grepl("^9606\\.", tmp$protein1))  # TRUE
# all(grepl("^9606\\.", tmp$protein2))  # TRUE
# # remove "9606." prefix
# tmp$protein1 <- gsub("^9606\\.", "", tmp$protein1)
# tmp$protein2 <- gsub("^9606\\.", "", tmp$protein2)
#
# # Map ENSP to HGNC symbols: use Dr. Steipe's mapping tool:
# load(file = file.path("inst", "extdata", "ensp2sym.RData"))
# tmp$protein1 <- ensp2sym[tmp$protein1]
# tmp$protein2 <- ensp2sym[tmp$protein2]
#
# # Validate initial mapping
# any(grepl("ENSP", tmp$protein1))  # Nope
# any(grepl("ENSP", tmp$protein2))  # None left here either
#
# # Clean duplicate edges (from Dr. Steipe)
# sPaste <- function(x, collapse = ":") {
#   return(paste(sort(x), collapse = collapse))
# }
# tmp$key <- apply(tmp[ , c("protein1", "protein2")], 1, sPaste) # takes a min
# length(tmp$key) # 35072
# length(unique(tmp$key)) # 17379
# tmp <- tmp[( ! duplicated(tmp$key)),
#            c("protein1", "protein2", "coexpression", "combined_score") ]
#
# # Remove NA nodes
# sum(is.na(tmp$protein1)) # 51
# sum(is.na(tmp$protein2)) # 153
# STRINGedges <- tmp[( ! is.na(tmp$protein1)) & ( ! is.na(tmp$protein2)), ] # 17175
#
# # Save the file
# save(STRINGedges, file = file.path(".", "data", "STRINGedges.RData"))
#
#
# #### PART 2 ####
# ##### Map the protein action dataset mappings ####
# tmp <- readr::read_tsv(file.path("./data", "9606.protein.actions.v11.0.txt"),
#                          skip = 1,
#                          col_names = c("protein1", "protein2", "mode", "action",
#                                         "is_directional", "a_is_acting", "combined_score")) # 11,759,454 rows
# # Keep "high confidence" interactions, and
# # remove "action" col since that information is duplicated in "mode" for our purposes
# tmp <- tmp[,c("protein1", "protein2", "mode",
#               "is_directional", "a_is_acting", "combined_score")]
# tmp <- tmp[tmp$combined_score >= 800, ]
#
# # remove "9606." prefix
# tmp$protein1 <- gsub("^9606\\.", "", tmp$protein1)
# tmp$protein2 <- gsub("^9606\\.", "", tmp$protein2)
#
# # Map ENSP to HGNC symbols: use Dr. Steipe's mapping tool:
# load(file = file.path("inst", "extdata", "ensp2sym.RData"))
# tmp$protein1 <- ensp2sym[tmp$protein1]
# tmp$protein2 <- ensp2sym[tmp$protein2]
#
# # Validate initial mapping
# any(grepl("ENSP", tmp$protein1))  # Nope
# any(grepl("ENSP", tmp$protein2))  # None left here either
#
# # Clean duplicate edges (from Dr. Steipe)
# sPaste <- function(x, collapse = ":") {
#   return(paste(sort(x), collapse = collapse))
# }
# tmp$key <- apply(tmp[ , c("protein1", "protein2", "mode")], 1, sPaste) # takes a min
# length(tmp$key) # 2031426
# length(unique(tmp$key)) # 548932
# tmp <- tmp[( ! duplicated(tmp$key)),
#            c("protein1", "protein2", "mode",
#              "is_directional", "a_is_acting", "combined_score") ]
#
# # Remove NA nodes
# sum(is.na(tmp$protein1)) # NUM
# sum(is.na(tmp$protein2)) # NUM
# STRINGactions <- tmp[( ! is.na(tmp$protein1)) & ( ! is.na(tmp$protein2)), ] # 545423
#
# # Save the file
# save(STRINGactions, file = file.path(".", "data", "STRINGactions.RData"))
#
# # How many of our high-confidence coexpression edges are in our protein actions dataset?
# validationEdges <- apply(STRINGedges[ , c("protein1", "protein2")], 1, sPaste)
# validationActions <- apply(STRINGactions[ , c("protein1", "protein2")], 1, sPaste)
# validate <- intersect(validationActions, validationEdges)
# View(validate)
# length(validate) # 6161 edges with logged protein actions
