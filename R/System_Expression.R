

exProf <- function(sym, hgnc = HGNC, ncol = 20) {
  # returns a set of numbers as a virtual expression profile, for
  # development purposes only.
  set.seed(which(hgnc$sym == sym))
  p <- as.vector(scale(runif(ncol)))
  set.seed(NULL)
  return(p)
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

System_Expression(sysname){
  myURL <- paste0("https://github.com/hyginn/",
                  "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
  load(url(myURL))  # loads HGNC data frame

  components <- fetchComponents(sysname)
  expressions <- sapply(components, exProf)
  meanExpression <- rowMeans(expressions)
  varExpression <- apply(expressions, 1, var)
  p.values <- pnorm(meanExpression, sd = sqrt(varExpression))
  p.values[p.values > 0.5] <- 1 - p.values[p.values > 0.5]
  UpConditions <- which(meanExpression > 0 & p.values < 0.5)
  DownConditions <- which(meanExpression < 0 & p.values < 0.5)

  SysExpression <- list(`Mean Expression` = meanExpression,
                        `Upregulated Conditions` = UpConditions,
                        `Downregulated Conditions` = DownConditions)

  return(SysExpression)
}
