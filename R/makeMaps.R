# makeMaps.R

#' Generate the genetic interpretation map of BioGrid GGI tags of system's physical interactions.
#'
#' @return The dataframe mapping BioGrid GGI tag to its interpretation assuming that the system's components also interact physically.
#' @examples
#' EMAP <- makeEMAP()
#'
#' @export
makeEMAP <- function() {
  geneticInteractions <- c(
    "Dosage Growth Defect",
    "Dosage Lethality",
    "Dosage Rescue",
    "Negative Genetic",
    "Phenotypic Enhancement",
    "Phenotypic Suppression",
    "Positive Genetic",
    "Synthetic Growth Defect",
    "Synthetic Haploinsufficiency",
    "Synthetic Lethality",
    "Synthetic Rescue"
  )

  effects <- c(
    "inhibited by",
    "inhibited by",
    "equivalent/activates",
    "cis-regulates",
    #negative genetic
    "trans-regulates",
    "cis-regulates",
    "trans-regulates",
    #positive genetic
    "equivalent/activates",
    "equivalent/activates",
    "equivalent/activates",
    "inhibited by"
  )

  notes <- c(
    "",
    "",
    "in a redundant pathway",
    "if protein 2 is inhibitor, protein 1 inhibits protein 2; protein 2 is an activator, protein 1 activates protein 2; potentially antagonistic",
    #negative genetic
    "if protein 2 is inhibitor, protein 1 activates protein 2; protein 2 is an activator, protein 1 inhibits protein 2",
    "if protein 2 is inhibitor, protein 1 inhibits protein 2; protein 2 is an activator, protein 1 activates protein 2",
    "potentially synergistic",
    #positive genetic
    "in a redundant pathway",
    "in a redundant pathway",
    "in a redundant pathway",
    ""
  )

  EMAP <- data.frame(geneticInt = geneticInteractions,
                     effect = effects,
                     notes = notes, stringsAsFactors = FALSE)

  return(EMAP)

}

#' Generate the genetic interpretation map of BioGrid GGI tags.
#'
#' @return The dataframe mapping BioGrid GGI tag to its interpretation.
#' @examples
#' GMAP <- makeGMAP()
#'
#' @export
makeGMAP <- function() {
  geneticInteractions <- c(
    "Dosage Growth Defect",
    "Dosage Lethality",
    "Dosage Rescue",
    "Negative Genetic",
    "Phenotypic Enhancement",
    "Phenotypic Suppression",
    "Positive Genetic",
    "Synthetic Growth Defect",
    "Synthetic Haploinsufficiency",
    "Synthetic Lethality",
    "Synthetic Rescue"
  )

  effects <- c(
    "negative-parallels",
    "negative-parallels",
    "positive-parallels",
    "synergizes with",
    #negative genetic
    "synergizes with",
    "antagonizes",
    "antagonizes",
    #positive genetic
    "synergizes with",
    "synergizes with",
    "synergizes with",
    "antagonizes"
  )

  GMAP <- data.frame(geneticInt = geneticInteractions,
                     effect = effects, stringsAsFactors = FALSE)

  return(GMAP)

}

# [END]
