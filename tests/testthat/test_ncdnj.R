library(EnvNJ)
context("NCD Approach")


## ----------------------------------------------------------- ##
#                 Testing the function ncd                      #
## ----------------------------------------------------------- ##
test_that("ncd() works properly", {

  skip_on_cran()

  a <- ncd(seq1 = "./data_t/Bos_taurus_prot.fasta",
           seq2 = "./data_t/Bison_bison_prot.fasta")
  b <- ncd(seq1 = "./data_t/Bos_taurus_dna.fasta",
           seq2 = "./data_t/Bison_bison_dna.fasta")

  expect_is(a, 'numeric')
  expect_equal(a, 0.0398, tolerance = 1e-5)
  expect_is(b, 'numeric')
  expect_equal(b, 0.2044, tolerance = 1e-5)
})

## ----------------------------------------------------------- ##
#               Testing the function ncdnj                      #
## ----------------------------------------------------------- ##
test_that("ncdnj() works properly", {

  skip_on_cran()

  a <- ncdnj(wd = "./data_t")

  expect_is(a, 'list')
  expect_is(a[[1]], 'matrix')
  expect_is(a[[2]], 'NULL')

})
