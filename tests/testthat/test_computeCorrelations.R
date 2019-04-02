library(testthat)
# Correlation graphs test (takes ~3min).
context("Correlation graphs test (takes ~170sec).")
expected <- list(
  "NLRIN" = list(recCount = 1225, corrSum = 145.111, goCorrSum = 674.092),
  "SLIGR" = list(recCount = 1830, corrSum = 147.04, goCorrSum = 1056.275),
  "PHALY" = list(recCount = 2145, corrSum = 101.0279, goCorrSum = 1348.659)
)
df <- BCB420.2019.ESA::computeCorrelations(TRUE)
expect_equal(expected["PHALY"][[1]], calcStats(df, "PHALY"), tolerance = 0.1)
expect_equal(expected["SLIGR"][[1]], calcStats(df, "SLIGR"), tolerance = 0.1)
expect_equal(expected["NLRIN"][[1]], calcStats(df, "NLRIN"), tolerance = 0.1)
