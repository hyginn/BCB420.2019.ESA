# test_makeSeq.R
# Author: Boris Steipe <https://orcid.org/0000-0002-1134-6758>
#
context("makeSeq")

# ==== BEGIN SETUP AND PREPARE =================================================
#
sampleSeq <- "ATGTATTTGTACCGGCGTTAA"

#
# ==== END SETUP AND PREPARE ===================================================

test_that("corrupt input generates errors",  {
  expect_error(makeSeq())
})

# CAUTION: this test sets the RNG seed.
test_that("a sample input produces the expected output",  {
  expect_equal(makeSeq(7, p = c(0.2, 0.4, 0.4, 0.2), seed = 112358),
               sampleSeq)
})


# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persistent construct that the test has created, except for
# stuff in tempdir().

# re-initialize RNG seed
set.seed(NULL)


#
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
