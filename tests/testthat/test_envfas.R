library(EnvNJ)
context("Env trees from fasta files")


## ----------------------------------------------------------- ##
#             Testing the function envfascpp                    #
## ----------------------------------------------------------- ##
test_that("envnj() works properly", {

  skip_on_cran()

  oldwd <- getwd()
  on.exit(oldwd)
  # setwd("/home/jca/Dropbox/Investigacion/Unpublished/EnvNJ/Env_trees/Local_Test_Linux")
  setwd("/Users/JCA/Dropbox/Investigacion/Unpublished/EnvNJ/Env_trees/Local_Test")
  envfascpp(path = ".", r = 10,
              exefile = ".",
              outfile = ".")

  expect_true(file.exists("v_Bison_bison_prot.txt"))
  expect_true(file.exists("v_Bos_taurus_prot.txt"))
  expect_true(file.exists("v_Dummy_species_prot.txt"))

  if (file.exists("v_Bison_bison_prot.txt")){
    file.remove("v_Bison_bison_prot.txt")
  }
  if (file.exists("v_Bos_taurus_prot.txt")){
    file.remove("v_Bos_taurus_prot.txt")
  }
  if (file.exists("v_Dummy_species_prot.txt")){
    file.remove("v_Dummy_species_prot.txt")
  }
  setwd(oldwd)
})

## ----------------------------------------------------------- ##
#             Testing the function vect2tree                    #
## ----------------------------------------------------------- ##
test_that("envnj() works properly", {

  skip_on_cran()

  a <- vect2tree(path = "./data_t")

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_is(a[[1]], 'matrix')
  expect_is(a[[2]], 'phylo')

})


