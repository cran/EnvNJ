library(EnvNJ)
context("Encoding OTUs")

## ---------------------------------------------- ##
#               Testing env.extract                #
## ---------------------------------------------- ##
test_that("env.extract() works properly",{

  skip_on_cran()
  skip_on_travis()
  seq1 <- "MFMINILMLIIPILLAVAFLTLVERKVLGYMQLRKGPNVVGPYGLLQPIADAIKLFIKEPLR"

  a <- env.extract(seq1, c = 11, r = 10)
  b <- env.extract(seq1, c = 3, r = 10)
  c <- env.extract(seq1, c = 60, r = 20)

  expect_is(a, 'character')
  expect_equal(nchar(a), 21)

  expect_is(b, 'character')
  expect_equal(nchar(b), 21)

  expect_is(c, 'character')
  expect_equal(nchar(c), 41)
})


## ---------------------------------------------- ##
#               Testing env.matrices               #
## ---------------------------------------------- ##
test_that("env.matrices() works properly",{

  a <- env.matrices(c('ANQRmCTPQ', 'LYPPmQTPC', 'XXGSmSGXX'))
  b <- env.matrices('ANQRmCTPQ')

  expect_is(a, 'list')
  expect_is(a[[1]], 'data.frame')
  expect_equal(nrow(a[[1]]), 2)
  expect_equal(ncol(a[[1]]), 9)
  expect_equal(as.character(a[[1]]$'0'), rep('m', 2))
  expect_is(a[[2]], 'data.frame')
  expect_equal(nrow(a[[2]]), 21)
  expect_equal(ncol(a[[2]]), 9)
  expect_equal(as.vector(a[[2]]$'0'), c(rep(0, 12), 2, rep(0, 8)))

  expect_is(b, 'list')
  expect_is(b[[1]], 'data.frame')
  expect_equal(nrow(b[[1]]), 1)
  expect_equal(ncol(b[[1]]), 9)
  expect_equal(as.character(b[[1]]$'0'), 'm')
  expect_is(b[[2]], 'data.frame')
  expect_equal(nrow(b[[2]]), 21)
  expect_equal(ncol(b[[2]]), 9)
  expect_equal(as.vector(b[[2]]$'0'), c(rep(0, 12), 1, rep(0, 8)))
})


## ----------------------------------------------------------- ##
#               Testing the function env.sp                     #
## ----------------------------------------------------------- ##
test_that("env.sp() works properly", {

  a <- env.sp(bovids, "Bos_taurus")
  b <- env.sp(bovids, "Syncerus_caffer", aa = c("M", "Y"))

  expect_is(a, 'list')
  expect_equal(length(a), 20)
  expect_equal(length(a[[1]]), 246)

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_equal(length(b[[1]]), 235)

})

## ----------------------------------------------------------- ##
#             Testing the function otu.vector                   #
## ----------------------------------------------------------- ##
test_that("otu.vector () works properly", {

  a <- otu.vector(envl = env.sp(bovids, "Bos_taurus", aa = "M"), aa = "M")
  b <- otu.vector(envl = env.sp(bovids, "Bison_bonasus", r = 3, aa = c("M", "Y")),
                  aa = c("M", "Y"))
  c <- otu.vector(envl = env.sp(bovids, "Bubalus_bubalis", r = 2))


  expect_is(a, 'matrix')
  expect_equal(dim(a), c(400, 1))
  expect_equal(sum(a), 4640)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(120, 2))
  expect_equal(sum(b[,1]), 1500)
  expect_equal(sum(b[,2]), 822)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(80, 20))
  expect_equal(sum(c[,1]), 1000)
  expect_equal(sum(c[,10]), 1304)
  expect_equal(sum(c[,20]), 752)
})

## ----------------------------------------------------------- ##
#             Testing the function otu.space                    #
## ----------------------------------------------------------- ##
test_that("otu.space() works properly", {

  a <- otu.space(bovids)

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(8000, 11))
  expect_equal(sum(a[,1]), 70600)
  expect_equal(sum(a[1,]), 191)
})


