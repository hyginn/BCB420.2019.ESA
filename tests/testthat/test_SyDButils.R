# test_SyDButils.R
# Author: Boris Steipe <https://orcid.org/0000-0002-1134-6758>
#
context("SyDButils")

# ==== BEGIN SETUP AND PREPARE =================================================
#

# Prepare test system DB

#
# ==== END SETUP AND PREPARE ===================================================

test_that("dummy",  {
  expect_equal(1,1)
})

# SyDBgetIDforKey()
# ToDo: test cases (+/- all)
#   - one val matches once in one column
#   - one val matches more than once in one column
#   - one val matches nowhere
#
#   - one val match once in multiple columns
#   - one val matches never in multiple columns
#   - one val match more than once in multiple columns
#
#   - multiple val all match once in one columns
#   - multiple val all match never in one column
#   - multiple val all match more than once in one column
#
#   - multiple val all match once in multiple columns
#   - multiple val all match never in multiple columns
#   - multiple val all match more than once in multiple columns

#   - some val match once in one column
#   - some val match more than once in one column
#   - some val match once in multiple columns

# SyDBgetIDforKey(all = FALSE, val = c("PHALY"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("nothg"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("concept"), att = c("moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("LC3-II"), att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("nothg"),  att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("concept"), att = c("name", "moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("autophagosome","LC3-II"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("foo","bar"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("concept", "RNA"), att = c("moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("autophagosome","LC3-II"), att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("foo","bar"), att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("Syntaxin 17", "RNA"), att = c("name", "moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("autophagosome","nothing"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("nothing", "RNA"), att = c("moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("autophagosome","nothing"), att = c("name", "code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = FALSE, val = c("Syntaxin 17", "nothing", "RNA"), att = c("name", "moleculeType"), tbl = "molecule", db = myDB)
#
# SyDBgetIDforKey(all = TRUE, val = c("PHALY"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("nothg"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("concept"), att = c("moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("LC3-II"), att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("nothg"),  att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("concept"), att = c("name", "moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("autophagosome","LC3-II"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("foo","bar"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("concept", "RNA"), att = c("moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("autophagosome","LC3-II"), att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("foo","bar"), att = c("code", "name"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("Syntaxin 17", "RNA"), att = c("name", "moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("autophagosome","nothing"), att = c("code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("nothing", "RNA"), att = c("moleculeType"), tbl = "molecule", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("autophagosome","nothing"), att = c("name", "code"), tbl = "system", db = myDB)
# SyDBgetIDforKey(all = TRUE, val = c("Syntaxin 17", "nothing", "RNA"), att = c("name", "moleculeType"), tbl = "molecule", db = myDB)


# SyDBTree
#
#  cat(SyDBTree(c("nonesuch"), myDB), sep = "\n")  # no entry
#  cat(SyDBTree(c("MARCH7"), myDB), sep = "\n")    # one leaf only
#  cat(SyDBTree(c("MARCH7", "FADD"), myDB), sep = "\n")    # two leaves
#  cat(SyDBTree(c("NF-kappa-B pathways"), myDB), sep = "\n")
#  cat(SyDBTree(c("NF-kappa-B", "NLRP3 up-regulation"), myDB), sep = "\n")
#  cat(SyDBTree(c("nonesuch", "NLRP3 up-regulation"), myDB), sep = "\n")
#  cat(SyDBTree(c("NF-kappa-B pathways", "NF-kappa-B"), myDB), sep = "\n")
#  cat(SyDBTree(c("PHALY"), myDB), sep = "\n")
#  cat(SyDBTree(c("SLIGR"), myDB), sep = "\n")
#  cat(SyDBTree(c("NLRIN"), myDB), sep = "\n")
#  cat(SyDBTree(c("HVGCR"), myDB), sep = "\n")
#  cat(SyDBTree(c("HVGCR"), myDB, MAX = 0), sep = "\n")
#  cat(SyDBTree(c("HVGCR"), myDB, MAX = 1), sep = "\n")
#  cat(SyDBTree(c("HVGCR"), myDB, MAX = 2), sep = "\n")
#  cat(SyDBTree(c("HVGCR"), myDB, MAX = 3), sep = "\n")



# ==== BEGIN TEARDOWN AND RESTORE ==============================================
# Remove every persistent construct that the test has created, except for
# stuff in tempdir().

# re-initialize RNG seed
set.seed(NULL)


#
# ==== END  TEARDOWN AND RESTORE ===============================================

# [END]
