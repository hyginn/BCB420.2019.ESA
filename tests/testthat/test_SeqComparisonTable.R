# test_SeqComparisonTable.R
# Author: Heewon (Judy) Lee <https://orcid.org/0000-0002-1134-6758>
# This is a validation test for SeqComparisonTable
#
context("SeqComparisonTable")

# ==== BEGIN SETUP AND PREPARE =================================================
#

# Choose a sample HGNC ID
hgnc <- "BECN1"
sampleSeqTable <- SeqComparisonTable(hgnc)
sampleSeqs <- c(sampleSeqTable[1,]$SequenceA, sampleSeqTable[1,]$SequenceB)

# Retrieve the sequence of each transcript from the PDB (retrieved from PDB directly)
# This is the sequence to the first chain of Beclin 1
BECN1201PdbSeq <- paste("DDSEQLQMELKELALEEERLIQELEDVEKNRKIVAENLEKVQAEAERLDQEE",
                  "AQYQREYSEFKRQQLELDDELKSVENQMRYAQTQLDKLKLE",
                  sep = "")

DYNLL2201PdbSeq <- paste("MSDRKAVIKNADMSEDMQQDAVDCATQAMEKYNIEKDIAAYIKKEFDKKYNP",
                        "TWHCIVGRNFGSYVTHETKHFIYFYLGQVAILLFKSG", sep = "")

compare <- c(BECN1201PdbSeq, DYNLL2201PdbSeq)

#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(SeqComparisonTable())
})

# CAUTION: this test sets the RNG seed.
test_that("a sample input produces the expected output",  {
  expect_equal(sampleSeqs, compare)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persistent construct that the test has created, except for
# stuff in tempdir().

# re-initialize RNG seed
set.seed(NULL)


#
# ==== END  TEARDOWN AND RESTORE ===============================================
#[END]
