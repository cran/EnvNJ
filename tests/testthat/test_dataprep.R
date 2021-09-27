library(EnvNJ)
context("Data Preparation")

## ----------------------------------------------------------- ##
#              Testing the function fastaconc                   #
## ----------------------------------------------------------- ##
test_that("fastaconc() works properly", {

  fastaconc(otus = c('A', 'B'), inputdir = "./data_t", out.file = "./HolaMundo.fasta" )

  expect_true(file.exists("./HolaMundo.fasta"))

  system("rm HolaMundo.fasta")
})


## ----------------------------------------------------------- ##
#               Testing the function df2fasta                   #
## ----------------------------------------------------------- ##
test_that("df2fasta() works properly", {

  skip_on_cran()

  df2fasta(df = bovids, out.file = "./xxx.fasta")

  expect_true(file.exists("./xxx.fasta"))

  system("rm xxx.fasta")

})

## ----------------------------------------------------------- ##
#               Testing the function d.phy2df                   #
## ----------------------------------------------------------- ##
test_that("d.phy2df() works properly", {

  skip_on_cran()

  a <- d.phy2df(phyfile = "./data_t/d_dummy.txt")
  b <- d.phy2df(phyfile = "./data_t/d_dummy.txt", as = 'data.frame')

  expect_is(a, 'matrix')
  expect_equal(dim(a), c(5,5))

  expect_is(b, 'data.frame')
  expect_equal(dim(a), c(5,5))
})

