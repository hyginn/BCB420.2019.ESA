# test_SyDButils.R

# Author: Denits Vasileva <https://orcid.org/0000-000?-????-????>

#

# ==== BEGIN SETUP AND PREPARE =================================================

#

library(testthat)

# Correlation graphs test (takes ~3min).

context("Correlation graphs test (takes ~170sec).")
# Prepare test system DB

expected <- list(
  "NLRIN" = list(recCount = 1225, corrSum = 145.111, goCorrSum = 674.092),
  "SLIGR" = list(recCount = 1830, corrSum = 147.04, goCorrSum = 1056.275),
  "PHALY" = list(recCount = 2145, corrSum = 101.0279, goCorrSum = 1348.659)
)



df <- computeCorrelations(TRUE)

calcPHALY <- BCB420.2019.ESA::calcStats(df, "PHALY")

calcSLIGR <- BCB420.2019.ESA::calcStats(df, "SLIGR")

calcNLRIN <- BCB420.2019.ESA::calcStats(df, "NLRIN")



tmp <- expected["PHALY"][[1]]

phalyExpRecCount <- as.numeric(tmp$recCount)

phalyExpCor <- as.numeric(tmp$corrSum)

phalyExpGo <- as.numeric(tmp$goCorrSum)



tmp <- expected["SLIGR"][[1]]

sligrExpRecCount <- as.numeric(tmp$recCount)

sligrExpCor <- as.numeric(tmp$corrSum)

sligrExpGo <- as.numeric(tmp$goCorrSum)



tmp <- expected["NLRIN"][[1]]

nlrinExpRecCount <- as.numeric(tmp$recCount)

nlrinExpCor <- as.numeric(tmp$corrSum)

nlrinExpGo <- as.numeric(tmp$goCorrSum)



phalyActRecCount <- as.numeric(calcPHALY$recCount)

phalyActCor <- as.numeric(calcPHALY$corrSum)

phalyActGo <- as.numeric(calcPHALY$goCorrSum)



sligrActRecCount <- as.numeric(calcSLIGR$recCount)

sligrActCor <- as.numeric(calcSLIGR$corrSum)

sligrActGo <- as.numeric(calcSLIGR$goCorrSum)



nlrinActRecCount <- as.numeric(calcNLRIN$recCount)

nlrinActCor <- as.numeric(calcNLRIN$corrSum)

nlrinActGo <- as.numeric(calcNLRIN$goCorrSum)



expect_equal(phalyExpRecCount, phalyActRecCount)

expect_equal(sligrExpRecCount,sligrActRecCount)

expect_equal(nlrinExpRecCount,nlrinActRecCount)



expect_equal(phalyExpCor,phalyActCor,tolerance = 0.01)

expect_equal(nlrinExpCor,nlrinActCor,tolerance = 0.01)

expect_equal(sligrExpCor,sligrActCor,tolerance = 0.01)



expect_equal(phalyExpGo,phalyActGo,tolerance = 0.01)

expect_equal(nlrinExpGo,nlrinActGo,tolerance = 0.01)

expect_equal(sligrExpGo,sligrActGo,tolerance = 0.01)

#

# ==== END SETUP AND PREPARE ===================================================



#==== BEGIN TEARDOWN AND RESTORE ==============================================

# Remove every persistent construct that the test has created, except for

# things in tempdir().

df <- NULL

#

# ==== END  TEARDOWN AND RESTORE ===============================================



# [END]

