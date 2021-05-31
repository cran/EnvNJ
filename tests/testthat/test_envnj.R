library(EnvNJ)
context("EnvNJ Trees")


## ----------------------------------------------------------- ##
#               Testing the function envnj                    #
## ----------------------------------------------------------- ##
test_that("envnj() works properly", {
  a <- envnj(bovids)
  b <- envnj(bovids, outgroup = "Pseudoryx_nghetinhensis")

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_is(a[[1]], 'matrix')
  expect_equal(dim(a[[1]]), c(11,11))

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_is(b[[1]], 'matrix')
  expect_equal(dim(b[[1]]), c(11,11))

})

## ----------------------------------------------------------- ##
#               Testing the function vcos                       #
## ----------------------------------------------------------- ##
test_that("vcos() works properly", {

  vectors1 <- list(a = c(0,0,1), b = c(1,1,0), c = c(0,0,1), d = c(1,1,1))
  a <- vcos(vectors1)
  vectors2 <- otu.space(bovids)
  b <- vcos(vectors2)

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(4,4))
  expect_equal(unname(a[1,2], force = TRUE), 0)
  expect_equal(unname(a[1,3], force = TRUE), 1)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(11,11))
  expect_equal(unname(b[1,2], force = TRUE), 0.995)
  expect_equal(unname(b[6,7], force = TRUE), 1)

})

## ----------------------------------------------------------- ##
#               Testing the function vdis                       #
## ----------------------------------------------------------- ##
test_that("vdis() works properly", {

  vectors <- list(a = c(0,0,1), b = c(1,1,0), c = c(0,0,1), d = c(1,1,1))
  a <- vdis(vcos(vectors))

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(4,4))
  expect_equal(unname(a[1,3], force = TRUE), 0)
})
