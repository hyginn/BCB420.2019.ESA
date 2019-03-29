# fetchComponents.R

#' \code{fetchComponents} get all gene symbols in a system.
#'
#' \code{fetchComponents} is deprecated. Use \code{SyDBgetSysSymbols()} instead.
#'
#' @param sys (character)  Name of a system
#'
#' @return (character) a vector of HGNC symbols
#' @examples
#' fetchComponents("PHALY")  # all PHALY system proteins
#' @export

fetchComponents <- function(sys) {
  # returns a fixed set of symbols.
  # Function stub for development purposes only.

  message(
    "Note: fetchComponents(...) is deprecated.",
    "Use SyDBgetSysSymbols(<database>, <system code>) instead.",
    sep = "\n"
  )

  if (sys == "PHALY") {
    if (sys == "PHALY") {
      s <-
        c(
          "AMBRA1",
          "ATG14",
          "ATP2A1",
          "ATP2A2",
          "ATP2A3",
          "BECN1",
          "BECN2",
          "BIRC6",
          "BLOC1S1",
          "BLOC1S2",
          "BORCS5",
          "BORCS6",
          "BORCS7",
          "BORCS8",
          "CACNA1A",
          "CALCOCO2",
          "CTTN",
          "DCTN1",
          "EPG5",
          "GABARAP",
          "GABARAPL1",
          "GABARAPL2",
          "HDAC6",
          "HSPB8",
          "INPP5E",
          "IRGM",
          "KXD1",
          "LAMP1",
          "LAMP2",
          "LAMP3",
          "LAMP5",
          "MAP1LC3A",
          "MAP1LC3B",
          "MAP1LC3C",
          "MGRN1",
          "MYO1C",
          "MYO6",
          "NAPA",
          "NSF",
          "OPTN",
          "OSBPL1A",
          "PI4K2A",
          "PIK3C3",
          "PLEKHM1",
          "PSEN1",
          "RAB20",
          "RAB21",
          "RAB29",
          "RAB34",
          "RAB39A",
          "RAB7A",
          "RAB7B",
          "RPTOR",
          "RUBCN",
          "RUBCNL",
          "SNAP29",
          "SNAP47",
          "SNAPIN",
          "SPG11",
          "STX17",
          "STX6",
          "SYT7",
          "TARDBP",
          "TFEB",
          "TGM2",
          "TIFA",
          "TMEM175",
          "TOM1",
          "TPCN1",
          "TPCN2",
          "TPPP",
          "TXNIP",
          "UVRAG",
          "VAMP3",
          "VAMP7",
          "VAMP8",
          "VAPA",
          "VPS11",
          "VPS16",
          "VPS18",
          "VPS33A",
          "VPS39",
          "VPS41",
          "VTI1B",
          "YKT6"
        )
    } else {
      s <- ""
    }

    return(s)
  }
}

# [END]
