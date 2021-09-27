## ---------- ngram.R -------------- ##
#                                     #
#    ngram                            #
#    ngraMatrix                       #
#    svdgram                          #
#                                     #
## --------------------------------- ##


## ----------------------------------------------------------- ##
#                   ngram(prot, k)                              #
## ----------------------------------------------------------- ##
#' Compute n-Gram Frequencies Vector
#' @description Computes the n-gram frequencies vector for a given protein.
#' @usage ngram(prot, k = 4)
#' @param prot a character string corresponding to the primary structure of the protein.
#' @param k a positive integer, between 1 and 5, indicating the k-mer of the words to be counted.
#' @details The one letter code for amino acids is used (capital).
#' @return A dataframe with two columns, the first one given the peptides and the second one the corresponding absolute frequency.
#' @references Stuart et al. Bioinformatics 2002; 18:100-108.
#' @seealso ngraMatrix(), ffp(), svdgram()
#' @examples ngram(bovids$Bos_taurus[1], k = 3)
#' @importFrom stringr str_count
#' @export

ngram <- function(prot, k = 4){

  aa <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
          'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  npe <- 20^k # number of peptides
  m <- matrix(rep(NA, npe * 2), nrow = npe)
  gram <- as.data.frame(m)
  names(gram) <- c('peptide', 'frequency')

  vector <- function(){
    v <- numeric(npe)
    for (j in 1:length(pept)){
      v[j] <- str_count(string = prot, pattern = paste("(?=", pept[j], ")", sep = ""))
    }
    return(v)
  }

  if (k == 1){

    gram$peptide <- pept <- aa
    gram$frequency <- vector()
    return(gram)

  } else if (k > 1){

    p <- aa # Base Case
    for (i in 1:(k-1)){
      p <- expand.grid(p, aa)
      p <- apply(p, 1, function(x) paste(x, collapse = ""))
    }

    gram$peptide <- pept <- p
    gram$frequency <- vector()
    return(gram)

  }
}


## ----------------------------------------------------------- ##
#              ngraMatrix(data, k, silent)                      #
## ----------------------------------------------------------- ##
#' Compute n-Gram Frequencies Dataframe
#' @description Computes the n-gram frequencies dataframe for the protein and species provides.
#' @usage ngraMatrix(data, k = 4, silent = FALSE)
#' @param data a dataframe with as many columns as species and one row per orthologous protein. The rows and columns must be named accordingly.
#' @param silent logical, set to FALSE to avoid loneliness.
#' @param k a positive integer, between 1 and 5, indicating the k-mer of the words to be counted.
#' @details The argument prot can be obtained using orth() and orth.seq().
#' @return A list with two dataframes. The first one with nsp * npr columns (nsp: number of species, npr: number of proteins per species) and npe rows (npe: number of peptides, 20 for n = 1, 400 for n = 2, 8000 for n = 3 and 160000 for n = 4). The entries of the dataframe are the number of times that the indicated peptide has been counted in the given protein. Orthologous proteins are in consecutive columns, thus the first nsp columns are the orthologous of protein 1 and so on. The second dataframe contains the Species Vector Sums (each vector describes one species).
#' @references Stuart et al. Bioinformatics 2002; 18:100-108.
#' @seealso ngram(), svdgram()
#' @examples ngraMatrix(bovids[,1:3], k = 2)
#' @export

ngraMatrix <- function(data, k = 4, silent = FALSE){

  aa <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
          'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
  nsp <- length(data) # number of species
  npr <- nrow(data) # number of proteins
  npe <- 20^k # number of peptides
  m <- matrix(rep(NA, npe * (nsp * npr)), nrow = npe)
  gram <- as.data.frame(m)
  names(gram)[1] <- 'peptide'

  p <- 0
  for (i in 1:npr){
    if (!silent){
      print(i)
    }
    for (j in 1:nsp){
      p <- p + 1
      t <- ngram(data[i,j], k = k)
      gram[, p+1] <- t$frequency
    }
  }
  gram$peptide <- t$peptide
  xx <- unlist(lapply(rownames(data), function(x) rep(x, nsp)))
  xx <- paste(xx, names(data), sep = "-")
  names(gram) <- c("peptide", xx)

  svs <- as.data.frame(matrix(rep(NA, npe * (nsp + 1)), nrow = npe))
  names(svs) <- c('peptide', names(data))
  svs$peptide <- gram$peptide
  counter <- 0
  for (i in 1:nsp){
    v <- numeric(nrow(svs))
    for (p in 1:npr){
      start <- 1 + nsp * (p - 1) + counter
      v <- v + gram[, start + 1]
    }
    svs[, i+1] <- v
    counter <- counter + 1
  }
  return(list(gram, svs))
}


## ----------------------------------------------------------- ##
#              svdgram(matrix, rank, species, SVS)              #
## ----------------------------------------------------------- ##
#' Compute Phylogenetic Trees Using an n-Gram and SVD Approach
#' @description Computes phylogenetic trees using an n-gram and SVD approach.
#' @usage svdgram(matrix, rank, species, SVS = TRUE)
#' @param matrix either a dataframe or a matrix where each row represents a property of a protein (for instance, the frequencies of tetrapeptides) and each column represents a different protein (or species).
#' @param rank a numeric array providing the ranks that want to be used to approach the data matrix using SVD.
#' @param species character array providing the species' names.
#' @param SVS logical. When the matrix passed as argument correspond to the peptide-protein matrix and SVS is set to TRUE, then the function will compute a matrix where the columns are the Species Vector Sums. Alternatively, if the matrix passed as argument is already a matrix where the columns encode for species, SVS should be set to FALSE.
#' @details When the matrix passed as argument is a matrix of peptide-protein, the function implement the method described by Stuart et al. 2002 (see references).
#' @return An object of class multiPhylo containing a tree for each rank value required.
#' @references Stuart et al. Bioinformatics 2002; 18:100-108.
#' @seealso ngraMatrix()
#' @examples
#' a <- ngraMatrix(bovids[, 1:4], k = 2)[[2]][, -1]
#' species <- names(a)
#' svdgram(matrix = a, rank = 4, species = species, SVS = FALSE)
#' @export

svdgram <- function(matrix, rank, species, SVS = TRUE){

  ## -- Check the selected ranks are allowed
  matrix <- as.matrix(matrix)
  sr <- sum(rank > min(dim(matrix)))
  if (sr > 0){
    stop("Rank should be below the matrix dimension!")
  }

  ## -- Check that the species names are provided
  if (length(species) == 0){
    stop("Please, provide species names")
  }

  trees <- vector(mode = "list", length(rank))

  singular <- svd(matrix)
  U <- singular$u
  s <- singular$d
  V <- singular$v

  Ak <- s[1] * (U[,1] %*% t(V[,1])) # best rank 1 matrix
  colnames(Ak) <- colnames(matrix)
  counter <- 0
  for (r in rank){
    counter <- counter + 1
    for (j in 2:r){
      print(paste("tree-rank: ", r, " ....", j, sep = ""))
      Ak <- Ak + (s[j] * (U[,j] %*% t(V[,j])))
    }
    if (SVS){
      Ak <- as.data.frame(Ak)
      Aksvs <- as.data.frame(matrix(rep(NA, dim(Ak)[1] * length(species)), ncol = length(species)))
      names(Aksvs) <- species
      for (i in 1:length(species)){
        sp <- species[i]
        at <- which(!is.na(stringr::str_extract(names(Ak), sp)))
        sub <- Ak[, at]
        Aksvs[i] <- apply(sub, 1, sum)
      }
      trees[[counter]] <- vtree(Aksvs)[[2]]
    } else {
      trees[[counter]] <- vtree(as.data.frame(Ak))[[2]]
    }
  }
  names(trees) <- paste("rank", rank, sep = "-")
  class(trees) <- "multiPhylo"
  return (trees)
}

