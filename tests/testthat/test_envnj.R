library(EnvNJ)
context("EnvNJ Trees")


## ----------------------------------------------------------- ##
#               Testing the function envnj                      #
## ----------------------------------------------------------- ##
test_that("envnj() works properly", {

  upgma <- envnj(bovids, clustering = 'upgma')
  a <- envnj(bovids, metric = 'euclidean')
  b <- envnj(bovids, outgroup = "Pseudoryx_nghetinhensis")
  c <- envnj(bovids, metric = 'manhattan')
  d <- envnj(bovids, metric = 'chebyshev')
  e <- envnj(bovids, metric = 'sorensen')
  f <- envnj(bovids, metric = 'soergel')
  g <- envnj(bovids, metric = 'lorentzian')
  h <- envnj(bovids, metric = 'kulczynski')
  i <- envnj(bovids, metric = 'canberra')
  j <- envnj(bovids, metric = 'non-intersection')
  k <- envnj(bovids, metric = 'wavehedges')
  l <- envnj(bovids, metric = 'czekanowski')
  m <- envnj(bovids, metric = 'motyka')
  n <- envnj(bovids, metric = 'cosine')
  o <- envnj(bovids, metric = 'jaccard')
  p <- envnj(bovids, metric = 'bhattacharyya')
  q <- envnj(bovids, metric = 'squared_chord')
  r <- envnj(bovids, metric = 'squared_chi')
  s <- envnj(bovids, metric = 'kullback-leibler')
  t <- envnj(bovids, metric = 'jeffreys')
  u <- envnj(bovids, metric = 'jensen-shannon')
  v <- envnj(bovids, metric = 'jensen_difference')
  w <- envnj(bovids, metric = 'hamming')
  x <- envnj(bovids, metric = 'mismatch')
  y <- envnj(bovids, metric = 'mismatchZero')
  z <- envnj(bovids, metric = 'binary')
  aa <- envnj(bovids, metric = 'taneja')
  bb <- envnj(bovids, metric = 'kumar-johnson')
  cc <- envnj(bovids, metric = 'avg')


  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_is(a[[1]], 'matrix')
  expect_equal(dim(a[[1]]), c(11,11))

  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_is(b[[1]], 'matrix')
  expect_equal(dim(b[[1]]), c(11,11))

  expect_is(c, 'list')
  expect_equal(length(c), 2)
  expect_is(c[[1]], 'matrix')
  expect_equal(dim(c[[1]]), c(11,11))

  expect_is(d, 'list')
  expect_equal(length(d), 2)
  expect_is(d[[1]], 'matrix')
  expect_equal(dim(d[[1]]), c(11,11))

  expect_is(e, 'list')
  expect_equal(length(e), 2)
  expect_is(e[[1]], 'matrix')
  expect_equal(dim(e[[1]]), c(11,11))

  expect_is(f, 'list')
  expect_equal(length(f), 2)
  expect_is(f[[1]], 'matrix')
  expect_equal(dim(f[[1]]), c(11,11))

  expect_is(g, 'list')
  expect_equal(length(g), 2)
  expect_is(g[[1]], 'matrix')
  expect_equal(dim(g[[1]]), c(11,11))

  expect_is(h, 'list')
  expect_equal(length(h), 2)
  expect_is(h[[1]], 'matrix')
  expect_equal(dim(h[[1]]), c(11,11))

  expect_is(i, 'list')
  expect_equal(length(i), 2)
  expect_is(i[[1]], 'matrix')
  expect_equal(dim(i[[1]]), c(11,11))

  expect_is(j, 'list')
  expect_equal(length(j), 2)
  expect_is(j[[1]], 'matrix')
  expect_equal(dim(j[[1]]), c(11,11))

  expect_is(k, 'list')
  expect_equal(length(k), 2)
  expect_is(k[[1]], 'matrix')
  expect_equal(dim(k[[1]]), c(11,11))

  expect_is(l, 'list')
  expect_equal(length(l), 2)
  expect_is(l[[1]], 'matrix')
  expect_equal(dim(l[[1]]), c(11,11))

  expect_is(m, 'list')
  expect_equal(length(m), 2)
  expect_is(m[[1]], 'matrix')
  expect_equal(dim(m[[1]]), c(11,11))

  expect_is(n, 'list')
  expect_equal(length(n), 2)
  expect_is(n[[1]], 'matrix')
  expect_equal(dim(n[[1]]), c(11,11))

  expect_is(o, 'list')
  expect_equal(length(o), 2)
  expect_is(o[[1]], 'matrix')
  expect_equal(dim(o[[1]]), c(11,11))

  expect_is(p, 'list')
  expect_equal(length(p), 2)
  expect_is(p[[1]], 'matrix')
  expect_equal(dim(p[[1]]), c(11,11))

  expect_is(q, 'list')
  expect_equal(length(q), 2)
  expect_is(q[[1]], 'matrix')
  expect_equal(dim(q[[1]]), c(11,11))

  expect_is(r, 'list')
  expect_equal(length(r), 2)
  expect_is(r[[1]], 'matrix')
  expect_equal(dim(r[[1]]), c(11,11))

  expect_is(s, 'list')
  expect_equal(length(s), 2)
  expect_is(s[[1]], 'matrix')
  expect_equal(dim(s[[1]]), c(11,11))

  expect_is(t, 'list')
  expect_equal(length(t), 2)
  expect_is(t[[1]], 'matrix')
  expect_equal(dim(t[[1]]), c(11,11))

  expect_is(u, 'list')
  expect_equal(length(u), 2)
  expect_is(u[[1]], 'matrix')
  expect_equal(dim(u[[1]]), c(11,11))

  expect_is(v, 'list')
  expect_equal(length(v), 2)
  expect_is(v[[1]], 'matrix')
  expect_equal(dim(v[[1]]), c(11,11))

  expect_is(w, 'list')
  expect_equal(length(w), 2)
  expect_is(w[[1]], 'matrix')
  expect_equal(dim(w[[1]]), c(11,11))

  expect_is(x, 'list')
  expect_equal(length(x), 2)
  expect_is(x[[1]], 'matrix')
  expect_equal(dim(x[[1]]), c(11,11))

  expect_is(y, 'list')
  expect_equal(length(y), 2)
  expect_is(y[[1]], 'matrix')
  expect_equal(dim(y[[1]]), c(11,11))

  expect_is(z, 'list')
  expect_equal(length(z), 2)
  expect_is(z[[1]], 'matrix')
  expect_equal(dim(z[[1]]), c(11,11))

  expect_is(aa, 'list')
  expect_equal(length(aa), 2)
  expect_is(aa[[1]], 'matrix')
  expect_equal(dim(aa[[1]]), c(11,11))

  expect_is(bb, 'list')
  expect_equal(length(bb), 2)
  expect_is(bb[[1]], 'matrix')
  expect_equal(dim(bb[[1]]), c(11,11))

  expect_is(cc, 'list')
  expect_equal(length(cc), 2)
  expect_is(cc[[1]], 'matrix')
  expect_equal(dim(cc[[1]]), c(11,11))

})


# ## ----------------------------------------------------------- ##
# #               Testing the function env.fasta                  #
# ## ----------------------------------------------------------- ##
# test_that("env.fasta() works properly", {
#
#   skip_on_cran()
#
#   a <- env.fasta(file = "./data_t/sample5.fasta")
#
#   expect_is(a, 'list')
#   expect_equal(length(a), 2)
#   expect_is(a[[1]], 'matrix')
#   expect_equal(dim(a[[1]]), c(5,5))
#   expect_is(a[[2]], 'phylo')
# })


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

  skip_on_cran()

  vectors <- list(a = c(0,0,1), b = c(1,1,0), c = c(0,0,1), d = c(1,1,1))
  a <- suppressWarnings(vdis(vcos(vectors)))

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(4,4))
  expect_equal(unname(a[1,3], force = TRUE), 0)
})
