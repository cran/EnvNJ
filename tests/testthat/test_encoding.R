library(EnvNJ)
context("Encoding OTUs")

## ----------------------------------------------- ##
#         Testing the function aa.at                #
## ----------------------------------------------- ##
test_that("aa.at() using named arguments with order", {

  skip_on_cran()
  skip_on_travis()

  a <- aa.at(at=28, target='P01009' , uniprot=TRUE)
  b <- aa.at(at=500, target='P01009' , uniprot=TRUE)

  expect_is(a, 'character')
  expect_equal(a, "Q")
  expect_is(b, 'NULL')

})

test_that("aa.at() using named arguments without order", {

  skip_on_cran()
  skip_on_travis()

  a <- aa.at(target='P01009' , uniprot=TRUE, at=28)
  b <- aa.at(at=500, uniprot=TRUE, target='P01009')


  expect_is(a, 'character')
  expect_equal(a, "Q")
  expect_equal(length(a), 1)
  expect_equal(nchar(a), 1)

  expect_is(b, 'NULL')

})

test_that("aa.at using uniprot = FALSE",{

  expect_equal(aa.at(6, "MARSSRAM", FALSE), 'R')
  expect_is(aa.at(10, "MARSSRAM", FALSE), 'NULL')

})


## ---------------------------------------------- ##
#               Testing env.extract                #
## ---------------------------------------------- ##
test_that("env.extract() works properly",{

  skip_on_cran()
  skip_on_travis()

  a <- env.extract('P01009', db = 'uniprot', c = 351+24, r = 10, ctr = 'random')
  b <- env.extract('P01009', db = 'uniprot', c = 351+24, r = 10, ctr = 'closest')
  c <- env.extract('P01009', db = 'uniprot', c = 7, r = 14)
  d <- env.extract('P01009', db = 'uniprot', c = 415, r = 20)
  e <- env.extract('P010091', db = 'uniprot', c = 15, r = 20)


  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equal(nchar(a$Positive), 21)
  expect_equal(nchar(a$Positive), nchar(a$Control))

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_equal(nchar(b$Positive), 21)
  expect_false(b$Positive == b$Control)

  expect_true(a$Positive == b$Positive)

  expect_is(c, 'list')
  expect_equal(length(c), 2)
  expect_true('X' %in% strsplit(c$Positive, split = "")[[1]])
  expect_equal(nchar(c$Control), 0)

  expect_is(d, 'list')
  expect_equal(length(d), 2)
  expect_true('X' %in% strsplit(d$Positive, split = "")[[1]])
  expect_equal(nchar(d$Control), 0)

  expect_is(e, 'NULL')
})


test_that("env.extract() works properly",{

  seq1 <- "ARSTVWXXWVTSRAYPILNMSSQQTTWWYYRTGFLIVSTHKRED"
  seq2 <- "ARSTVWXXWVTSRAYPILNMSSQQTTWWYYMRTGFLIVSTHKRED"
  a <- env.extract(seq1,  c = 20, r = 5, ctr = 'random')
  b <- env.extract(seq2, c = 20, r = 5, ctr = 'random')
  c <- env.extract(seq2,  c = 20, r = 5, ctr = 'closest')
  d <- env.extract(seq2,  c = 20, r = 5, ctr = 'random', exclude = 31)

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equal(nchar(a$Positive), 11)
  expect_gt(nchar(a$Positive), nchar(a$Control))

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_equal(nchar(b$Positive), 11)
  expect_false(b$Positive == b$Control)
  expect_true(a$Positive == b$Positive)

  expect_is(c, 'list')
  expect_equal(length(c), 2)
  expect_true('W' %in% strsplit(c$Control, split = "")[[1]])
  expect_equal(nchar(c$Positive), nchar(c$Control))

  expect_is(d, 'list')
  expect_equal(length(d), 2)
  expect_true('S' %in% strsplit(d$Positive, split = "")[[1]])
  expect_equal(nchar(d$Control), 0)
})


## ---------------------------------------------- ##
#               Testing env.matrices               #
## ---------------------------------------------- ##
test_that("env.matrices() works properly",{

  a <- env.matrices(c('ANQRmCTPQ', 'LYPPmQTPC', 'XXGSmSGXX'))

  expect_is(a, 'list')
  expect_is(a[[1]], 'data.frame')
  expect_equal(nrow(a[[1]]), 3)
  expect_equal(ncol(a[[1]]), 9)
  expect_equal(as.character(a[[1]]$'0'), rep('m', 3))
  expect_is(a[[2]], 'data.frame')
  expect_equal(nrow(a[[2]]), 21)
  expect_equal(ncol(a[[2]]), 9)
  expect_equal(as.vector(a[[2]]$'0'), c(rep(0, 12), 3, rep(0, 8)))
})


## ----------------------------------------------------------- ##
#               Testing the function env.sp                     #
## ----------------------------------------------------------- ##
test_that("env.sp() works properly", {

  a <- env.sp(bovids, "Bos_taurus")
  b <- env.sp(bovids, "Bos_taurus", remove.init = FALSE)
  c <- env.sp(bovids, "Bos_taurus", r = 3)
  d <- env.sp(bovids, "Syncerus_caffer")
  e <- env.sp(bovids, "Bos_taurus", aa = c("M", "Y"))

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 3517)
  expect_equal(ncol(a), 6)
  expect_equal(nchar(a$env[1]), 21)
  expect_equal(a$species[1], "Bos_taurus")

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 3530)
  expect_equal(ncol(b), 6)
  expect_equal(nchar(b$env[1]), 21)
  expect_equal(b$species[1], "Bos_taurus")

  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 3699)
  expect_equal(ncol(c), 6)
  expect_equal(nchar(c$env[1]), 7)
  expect_equal(c$species[1], "Bos_taurus")

  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 3517)
  expect_equal(ncol(d), 6)
  expect_equal(nchar(d$env[1]), 21)
  expect_equal(d$species[1], "Syncerus_caffer")

  expect_is(e, 'data.frame')
  expect_equal(nrow(e), 356)
  expect_equal(ncol(e), 6)
  expect_equal(nchar(e$env[1]), 21)
  expect_equal(e$species[1], "Bos_taurus")
})

## ----------------------------------------------------------- ##
#             Testing the function otu.vector                   #
## ----------------------------------------------------------- ##
test_that("otu.vector () works properly", {

  a <- otu.vector(df = env.sp(bovids, "Bos_taurus"), aa = "M")
  b <- otu.vector(df = env.sp(bovids, "Bison_bonasus", r = 3), aa = c("M", "Y"))
  c <- otu.vector(df = env.sp(bovids, "Bubalus_bubalis", r = 2))

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(400, 1))
  expect_equal(sum(a), 4620)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(120, 2))
  expect_equal(sum(b[,1]), 1494)
  expect_equal(sum(b[,2]), 810)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(80, 20))
  expect_equal(sum(c[,1]), 1000)
  expect_equal(sum(c[,10]), 1300)
  expect_equal(sum(c[,20]), 752)

})

## ----------------------------------------------------------- ##
#             Testing the function otu.space                    #
## ----------------------------------------------------------- ##
test_that("otu.space() works properly", {

  a <- otu.space(bovids)

  expect_is(a, 'data.frame')
  expect_equal(dim(a), c(8000, 11))

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
