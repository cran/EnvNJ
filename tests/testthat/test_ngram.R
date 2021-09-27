library(EnvNJ)
context("n-Gram")


## ----------------------------------------------------------- ##
#               Testing the function ngram                      #
## ----------------------------------------------------------- ##
test_that("ngram() works properly", {

  s <- "MSAGARRRPR"
  a <- ngram(s, k = 1)
  b <- ngram(s, k = 2)
  c <- ngram(s, k = 3)
  d <- ngram(s, k = 4)

  expect_is(a, 'data.frame')
  expect_equal(length(a), 2)
  expect_equal(nrow(a), 20)
  expect_equal(sum(a$frequency), 10)

  expect_is(b, 'data.frame')
  expect_equal(length(b), 2)
  expect_equal(nrow(b), 400)
  expect_equal(sum(b$frequency), 9)

  expect_is(c, 'data.frame')
  expect_equal(length(c), 2)
  expect_equal(nrow(c), 8000)
  expect_equal(sum(c$frequency), 8)

  expect_is(d, 'data.frame')
  expect_equal(length(d), 2)
  expect_equal(nrow(d), 160000)
  expect_equal(sum(d$frequency), 7)

})


## ----------------------------------------------------------- ##
#           Testing the function ngraMatrix                     #
## ----------------------------------------------------------- ##
test_that("ngraMatrix() works properly", {

  prot <- data.frame(sp1 = c("MSPGASRGPR", "MKMADRSGKI", "MACTIQKAEA"),
                     sp2 = c("MSAGARRRPR", "MADRGKLVPG", "MCCAIQKAEA"))
  a <- ngraMatrix(prot, k = 1)
  a1 <- a[[1]]
  a2 <- a[[2]]

  expect_is(a, 'list')
  expect_is(a1, 'data.frame')
  expect_is(a2, 'data.frame')
  expect_equal(nrow(a1), 20)
  expect_equal(length(a1), 7)
  expect_equal(nrow(a2), 20)
  expect_equal(length(a2), 3)
  expect_equal(sum(apply(a1[-1], 2, sum) == rep(10,6)), 6)
  expect_equal(sum(apply(a2[-1], 2, sum) == rep(30,2)), 2)

  b <- ngraMatrix(prot, k = 2)
  b1 <- b[[1]]
  b2 <- b[[2]]

  expect_is(b, 'list')
  expect_is(b1, 'data.frame')
  expect_is(b2, 'data.frame')
  expect_equal(nrow(b1), 400)
  expect_equal(length(b1), 7)
  expect_equal(nrow(b2), 400)
  expect_equal(length(b2), 3)
  expect_equal(sum(apply(b1[-1], 2, sum) == rep(9,6)), 6)
  expect_equal(sum(apply(b2[-1], 2, sum) == rep(27,2)), 2)

  c <- ngraMatrix(prot, k = 3)
  c1 <- c[[1]]
  c2 <- c[[2]]

  expect_is(c, 'list')
  expect_is(c1, 'data.frame')
  expect_is(c2, 'data.frame')
  expect_equal(nrow(c1), 8000)
  expect_equal(length(c1), 7)
  expect_equal(nrow(c2), 8000)
  expect_equal(length(c2), 3)
  expect_equal(sum(apply(c1[-1], 2, sum) == rep(8,6)), 6)
  expect_equal(sum(apply(c2[-1], 2, sum) == rep(24,2)), 2)

  d <- ngraMatrix(prot, k = 4)
  d1 <- d[[1]]
  d2 <- d[[2]]

  expect_is(d, 'list')
  expect_is(d1, 'data.frame')
  expect_is(d2, 'data.frame')
  expect_equal(nrow(d1), 160000)
  expect_equal(length(d1), 7)
  expect_equal(nrow(d2), 160000)
  expect_equal(length(d2), 3)
  expect_equal(sum(apply(d1[-1], 2, sum) == rep(7,6)), 6)
  expect_equal(sum(apply(d2[-1], 2, sum) == rep(21,2)), 2)
})


## ----------------------------------------------------------- ##
#                Testing the function svdgram                   #
## ----------------------------------------------------------- ##
test_that("svdgram() works properly", {

  a <- ngraMatrix(bovids, k = 2)[[1]][, -1]
  b <- ngraMatrix(bovids, k = 2)[[2]][, -1]
  species <- names(b)

  ta <- svdgram(matrix = a, rank = c(50, 143), species = species, SVS = TRUE)
  tb <- svdgram(matrix = b, rank = 11, species = species, SVS = FALSE)

  expect_is(ta, "multiPhylo")
  expect_equal(ta[[1]]$Nnode, 9)
  expect_equal(sum(ta[[1]]$tip.label == species), 11)
  expect_equal(names(ta[1]), "rank-50")

  expect_is(tb, "multiPhylo")
  expect_equal(tb[[1]]$Nnode, 9)
  expect_equal(sum(tb[[1]]$tip.label == species), 11)
  expect_equal(names(tb[1]), "rank-11")
})

