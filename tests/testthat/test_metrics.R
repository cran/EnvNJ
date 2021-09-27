library(EnvNJ)
context("Pairwise Vector Dissimilarities")

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
#               Testing the function cos2dis                    #
## ----------------------------------------------------------- ##
test_that("cos2dis() works properly", {

  vectors <- list(a = c(0,0,1), b = c(1,1,0), c = c(0,0,1), d = c(1,1,1))
  a <- cos2dis(vcos(vectors))

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(4,4))
  expect_equal(unname(a[1,3], force = TRUE), 0)
})


## ----------------------------------------------------------- ##
#               Testing the function cos2dis                    #
## ----------------------------------------------------------- ##
test_that("metrics() works properly", {

  ## ------------------------------------------------ L_p family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'euclidean')
  b <- metrics(vset, 'chebyshev')
  c <- metrics(vset, 'manhattan')
  d <- metrics(vset, 'minkowski', p = 3)

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equal(a[1,2], 17.52142, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equal(b[1,2], 16)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(3,3))
  expect_equal(c, t(c))
  expect_equal(c[1,2], 29)

  expect_is(d, 'matrix')
  expect_equal(dim(d), c(3,3))
  expect_equal(d, t(d))
  expect_equal(d[1,2], 16.277704)

})

test_that("metrics() works properly", {
  ## ------------------------------------------------ L_1 family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'sorensen')
  b <- metrics(vset, 'soergel')
  c <- metrics(vset, 'kulczynski')
  d <- metrics(vset, 'lorentzian')
  e <- metrics(vset, 'canberra')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equal(a[1,2], 0.5471698, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equal(b[1,2], 0.7073171, tolerance = 1e-5)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(3,3))
  expect_equal(c, t(c))
  expect_equal(c[1,2], 2.416659, tolerance = 1e-5)

  expect_is(d, 'matrix')
  expect_equal(dim(d), c(3,3))
  expect_equal(d, t(d))
  expect_equal(d[1,2], 8.313852, tolerance = 1e-5)

  expect_is(e, 'matrix')
  expect_equal(dim(e), c(3,3))
  expect_equal(e, t(e))
  expect_equal(e[1,2], 4.533333, tolerance = 1e-5)
})

test_that("metrics() works properly", {
  ## --------------------------------------- Intersection family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'non-intersection')
  b <- metrics(vset, 'wavehedges')
  c <- metrics(vset, 'czekanowski')
  d <- metrics(vset, 'motyka')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equal(a[1,2], 0.3365385, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equal(b[1,2], 4.695652, tolerance = 1e-5)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(3,3))
  expect_equal(c, t(c))
  expect_equal(c[1,2], 0.5471698)

  expect_is(d, 'matrix')
  expect_equal(dim(d), c(3,3))
  expect_equal(d, t(d))
  expect_equal(d[1,2], 0.7735849)
})

test_that("metrics() works properly", {
  ## ------------------------------------- Inner product family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'cosine')
  b <- metrics(vset, 'jaccard')


  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equivalent(a[1,2], 0.05604139, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equal(b[1,2], 0.6329897, tolerance = 1e-5)
})

test_that("metrics() works properly", {
  ## ------------------------------------- Squared-chord family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'bhattacharyya')
  b <- metrics(vset, 'squared_chord')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equivalent(a[1,2], 0.2539953, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equal(b[1,2], 17.62284, tolerance = 1e-5)
})

test_that("metrics() works properly", {
  ## ------------------------------------- Squared Chi family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'squared_chi')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equivalent(a[1,2], 21.53333, tolerance = 1e-5)
})


test_that("metrics() works properly", {
  ## ---------------------------------- Shannon entropy family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'kullback-leibler')
  b <- metrics(vset, 'jeffreys')
  c <- metrics(vset, 'jensen-shannon')
  d <- metrics(vset, 'jensen_difference')


  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equivalent(a[1,2], 2.666679, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equivalent(b[1,2], 3.949217, tolerance = 1e-5)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(3,3))
  expect_equal(c, t(c))
  expect_equivalent(c[1,2], 0.1655989, tolerance = 1e-5)

  expect_is(d, 'matrix')
  expect_equal(dim(d), c(3,3))
  expect_equal(d, t(d))
  expect_equivalent(d[1,2], 0.1655989, tolerance = 1e-5)
})


test_that("metrics() works properly", {
  ## ----------------------------------------- Mismatch family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'hamming')
  b <- metrics(vset, 'mismatch')
  c <- metrics(vset, 'mismatchZero')
  d <- metrics(vset, 'binary')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equivalent(a[1,2], 0.7142857, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equivalent(b[1,2], 5, tolerance = 1e-5)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(3,3))
  expect_equal(c, t(c))
  expect_equivalent(c[1,2], 0.7142857, tolerance = 1e-5)

  expect_is(d, 'matrix')
  expect_equal(dim(d), c(3,3))
  expect_equal(d, t(d))
  expect_equivalent(d[1,2], 0.5714286, tolerance = 1e-5)
})


test_that("metrics() works properly", {
  ## -------------------------------------- Combinations family
  vset <- matrix(c(4, 23,  4,  0,  5,  1,  3,
                   0,  7,  4,  1,  0,  1,  0,
                   0,  0,  1,  1,  0,  0,  1), ncol = 3, byrow = FALSE)
  a <- metrics(vset, 'taneja')
  b <- metrics(vset, 'kumar-johnson')
  c <- metrics(vset, 'avg')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(3,3))
  expect_equal(a, t(a))
  expect_equivalent(a[1,2], 1.767442, tolerance = 1e-5)

  expect_is(b, 'matrix')
  expect_equal(dim(b), c(3,3))
  expect_equal(b, t(b))
  expect_equivalent(b[1,2], 41.91446, tolerance = 1e-5)

  expect_is(c, 'matrix')
  expect_equal(dim(c), c(3,3))
  expect_equal(c, t(c))
  expect_equivalent(c[1,2], 22.5, tolerance = 1e-5)
})
