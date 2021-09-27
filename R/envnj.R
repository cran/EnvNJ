## ----------- envnj.R ------------- ##
#                                     #
#    envnj                            #
#    env.fasta                        #
#    vdis                             #
#                                     #
## --------------------------------- ##


## ----------------------------------------------------------- ##
#        envnj(data, r, aa, metric, clustering, outgroup)       #
## ----------------------------------------------------------- ##
#' Build Trees Based on the Environment Around the Indicated Amino Acid(s)
#' @description Builds trees based on the environment around the indicated amino acid(s).
#' @usage envnj(data, r = 10, aa = 'all', metric = "cosine", clustering = "nj", outgroup = 'any')
#' @param data input data must be a dataframe where each row corresponds to a protein sequence and each column to a species.
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) to be used to encoded the species.
#' @param metric character string indicating the metric (see metrics() to see the methods allowed).
#' @param clustering string indicating the clustering method, either "nj" or "upgma".
#' @param outgroup when a rooted tree is desired, it indicates the species to be used as outgroup.
#' @details This function builds alignment-independent phylogenetic trees.
#' @return A list with two objects, the first one is an inter-species distance matrix. The second one is an object of class 'phylo'.
#' @seealso otu.space(), metrics()
#' @examples \donttest{
#' data(bovids)
#' envnj(bovids[, 7:11], aa = "all", outgroup = "Pseudoryx_nghetinhensis")
#' }
#' @importFrom ape nj
#' @importFrom ape root
#' @importFrom ape plot.phylo
#' @importFrom phangorn upgma
#' @export

envnj <- function(data, r = 10, aa = "all", metric = "cosine", clustering = "nj", outgroup = 'any'){

  space <- otu.space(data = data, r = r, aa = aa)
  d <- metrics(vset = space, method = metric)
  if (clustering == 'nj'){
    t <- ape::nj(d)
  } else if (clustering == 'upgma'){
    t <- phangorn::upgma(d)
  } else {
    warning("NJ method has been selected")
  }

  if (outgroup[1] != "any"){
    t <- ape::root(t, outgroup = outgroup)
  }
  ape::plot.phylo(t, use.edge.length = FALSE, cex = 0.5)
  return(list(d, t))
}

## ---------------------------------------------------- ##
##                    env.fasta()                       ##
## ---------------------------------------------------- ##
#' Build Trees Based on the Environment Around the Indicated Amino Acid(s)
#' @description Builds trees based on the environment around the indicated amino acid(s).
#' @usage env.fasta(file, r = 10, aa = 'all', out.file = 'any')
#' @param file path to the single multispecies fasta file to be used as input.
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) to be used to encoded the species.
#' @param out.file path and name of output file. Only if intermediate results data want to be saved (see details).
#' @details This function builds alignment-independent phylogenetic trees. The input data is a fasta file. When an out.file path is provided, the environment sequences of each species and the vector representing each species are saved in the path provided.
#' @return A list with two objects, the first one is an inter-species distance matrix. The second one is an object of class 'phylo'.
#' @seealso envnj(), fastaconc()
#' @examples \dontrun{env.fasta(file = "./data_t/sample5.fasta")}
#' @importFrom ape nj
#' @importFrom ape plot.phylo
#' @importFrom seqinr read.fasta
#' @export

env.fasta <- function(file, r = 10, aa = "all", out.file = 'any'){

  if (out.file != 'any'){
    dir.create(path = out.file)
    dir.create(paste(out.file, "/Extracted_Environments_r", r, sep = ""))
    dir.create(paste(out.file, "/Species_Vectors_r", r, sep = ""))
  }

  if (aa == "all"){
    naa <- 20
  } else {
    naa <- length(aa)
  }
  # Read the input fasta file:
  f <- seqinr::read.fasta(file, seqtype = 'AA', as.string = TRUE)
  species <- names(f)

  space <- matrix(rep(NA, naa * 40*r * length(species)), ncol = length(species))
  colnames(space) <- species

  for (i in 1:length(species)){
    print(paste(i, "  -----  ", species[i]))

    pdata <- data.frame(sp = as.character(f[[i]]))
    names(pdata) <- species[i]

    env <- EnvNJ::env.sp(data = pdata, sp = species[i], r = r, silent = FALSE)
    v <- EnvNJ::otu.vector(env)
    v <- as.vector(v)

    if (out.file != "any"){
      suppressWarnings(
        save(env, file = paste(out.file, "/Extracted_Environments_r", r, "/env_",
                               species[i], ".Rda", sep = ""))
      )
      suppressWarnings(
        save(v, file = paste(out.file, "/Species_Vectors_r", r, "/v_",
                             species[i], ".Rda", sep = ""))
      )
    }

    space[, i] <- v
  }
  space <- as.data.frame(space)

  if (out.file != 'any'){
    save(space,
         file = paste(out.file, "/Species_Vector_Subspace_r", r, ".Rda", sep = ""))

    print(paste("Results saved at ", out.file, sep = ""))
  }

  cosDist <- EnvNJ::vcos(space, silent = TRUE, digits = 6)
  d <- cos2dis(cosDist)
  d[is.na(d)] <- 0
  d <- d + t(d)
  t <- ape::nj(d)

  ape::plot.phylo(t, use.edge.length = FALSE, cex = 0.5)

  return(list(d, t))
}

## ---------------------------------------------------- ##
##                    vdis()                            ##
## ---------------------------------------------------- ##
#' Compute Pairwise Distances Between Vectors
#' @description Computes pairwise distances between vectors.
#' @usage vdis(cos)
#' @param cos a square upper triangular matrix where cos(i,j) is the cosine between the vector i and j.
#' @details Cosines are standard measure of vector similarity, and can be converted into distance by dij = -log( (1 + cos(i,j) )/2).
#' @return A triangular matrix with the distances.
#' @seealso vcos()
#' @examples
#' data(bovids)
#' vectors = otu.space(bovids[, 7:11])
#' cosData = vcos(vectors)
#' disData = suppressWarnings(vdis(cosData))
#' @export


vdis <- function(cos){

  .Deprecated("cos2dis")
  # Check that cos is a square matrix
  if (!is.matrix(cos)){
    stop("The input must be a matrix")
  } else if (length(unique(dim(cos))) != 1){
    stop("The input matrix must be square")
  }

  d <- cos
  n <- dim(d)[1]
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      d[i,j] <- -log((1 + cos[i,j])/2)
    }
  }
  return(d)
}
