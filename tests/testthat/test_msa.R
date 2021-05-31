library(EnvNJ)
context("MSA Related")

## ----------------------------------------------------------- ##
#             Testing the function msa.merge                    #
## ----------------------------------------------------------- ##
test_that("msa.merge() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- msa.merge(bovids)

  expect_is(a, 'data.frame')
  expect_equal(dim(a), c(11, 3791))
  expect_true("Bos_taurus" %in% rownames(a))

})

## ----------------------------------------------------------- ##
#               Testing the function msa.tree                   #
## ----------------------------------------------------------- ##
test_that("msa.tree() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- msa.tree(bovids)

  expect_is(a, "list")
  expect_equal(length(a), 3)
  expect_is(a[[1]], 'data.frame')
  expect_equal(dim(a[[1]]), c(11, 3791))
  expect_is(a[[2]], 'matrix')
  expect_equal(dim(a[[2]]), c(11,11))
  expect_is(a[[3]], 'phylo')
})

