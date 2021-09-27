library(EnvNJ)
context("Ancillary tools")

## ----------------------------------------------- ##
#         Testing the function aa.at                #
## ----------------------------------------------- ##
test_that("aa.at() works properly", {

  a <- aa.at(at = 28, target= bovids[1,1])

  expect_is(a, 'character')
  expect_equal(a, "L")

})


## ---------------------------------------------- ##
#        Testing the function aa.comp              #
## ---------------------------------------------- ##
test_that("the function aa.comp() works properly", {

  skip_on_cran()
  skip_on_travis()

  a <- aa.comp("ACEEAGRKDNW", uniprot = FALSE)
  b <- aa.comp("P01009")
  c <- aa.comp('P010091')

  expect_is(a, 'data.frame')
  expect_equal(dim(a), c(20, 2))

  if (!is.null(b)){
    expect_is(b, 'data.frame')
    expect_equal(dim(b), c(20, 2))
    expect_equal(sum(b$frequency), 418)
  }

  expect_is(c, 'NULL')
})


## ----------------------------------------------------------- ##
#               Testing the function aaf                        #
## ----------------------------------------------------------- ##
test_that("aaf() works properly", {

  a <- aaf(bovids)

  expect_is(a, 'data.frame')
  expect_equal(dim(a), c(20, 11))
  expect_equal(sum(a[,1]), 3790)

})
