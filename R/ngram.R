## ---------- ngram.R -------------- ##
#                                     #
#    aa.comp                          #
#    ngram                            #
#    ngraMatrix                       #
#    svdgram                          #
#    vtree                            #
#                                     #
## --------------------------------- ##


## --------------------------------------------------------------- ##
#             aa.comp(target, uniprot = TRUE)                       #
## --------------------------------------------------------------- ##
#' Amino Acid Composition
#' @description Returns a table with the amino acid composition of the target protein.
#' @usage aa.comp(target, uniprot = TRUE)
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param uniprot logical, if TRUE the argument 'target' should be an ID.
#' @return Returns a dataframe with the absolute frequency of each type of residue found in the target peptide.
#' @examples aa.comp('MPSSVSWGILLLAGLCCLVPVSLAEDPQGDAAQK', uniprot = FALSE)
#' @importFrom bio3d get.seq
#' @export

aa.comp <- function(target, uniprot = TRUE){

  if (uniprot == TRUE){
    seq <- tryCatch(
      {
        bio3d::get.seq(id = target)$ali
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(seq)){
      message("Sorry, get.seq failed")
      return(NULL)
    }

    seq <- paste(seq, collapse = "")
    if (file.exists("seqs.fasta")){
      file.remove("seqs.fasta")
    }
    id <- target

  } else {
    seq <- target
    id <- 'user-provided sequence'
  }

  aai <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
           "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  output <- data.frame(aa = aai, frequency = NA)

  for (aa in output$aa){
    t <- gregexpr(aa, seq)[[1]]
    if (t[1] == -1){
      output$frequency[which(output$aa == aa)] <- 0
    } else {
      output$frequency[which(output$aa == aa)] <- length(t)
    }
  }

  attr(output, 'seq') <- target
  return(output)
}


## ----------------------------------------------------------- ##
#                   ngram(prot, n)                              #
## ----------------------------------------------------------- ##
#' Compute n-Gram Frequencies Vector
#' @description Computes the n-gram frequencies vector for a given protein.
#' @usage ngram(prot, n = 1)
#' @param prot a character string corresponding to the primary structure of the protein.
#' @param n a positive integer between 1 and 4.
#' @details The one letter code for amino acids is used (capital).
#' @return A dataframe with two columns, the first one given the peptides and the second one the corresponding absolute frequency.
#' @references Stuart et al. Bioinformatics 2002; 18:100-108.
#' @seealso ngraMatrix(), svdgram()
#' @examples ngram(bovids$Bos_taurus[1], n = 3)
#' @importFrom stringr str_count
#' @export

ngram <- function(prot, n = 1){

  aa <- c("A","C","D","E","F","G","H","I","K","L","M",
          "N","P","Q","R","S","T","V","W","Y")
  npe <- 20^n # number of peptides
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

  if (n == 1){

    gram$peptide <- aa
    gram$frequency <- aa.comp(prot, uniprot = FALSE)$frequency
    return(gram)

  } else if (n == 2){

    pept <- expand.grid(aa, aa)
    pept <- apply(pept, 1, function(x) paste(x, collapse = ""))
    gram$peptide <- pept
    gram$frequency <- vector()
    return(gram)

  } else if  (n == 3){

    p <- expand.grid(aa, aa)
    p <- apply(p, 1, function(x) paste(x, collapse = ""))
    pept <- expand.grid(p, aa)
    pept <- apply(pept, 1, function(x) paste(x, collapse = ""))
    gram$peptide <- pept
    gram$frequency <- vector()
    return(gram)

  } else if (n == 4){

    p2 <- expand.grid(aa, aa)
    p2 <- apply(p2, 1, function(x) paste(x, collapse = ""))
    p3 <- expand.grid(p2, aa)
    p3 <- apply(p3, 1, function(x) paste(x, collapse = ""))
    p4 <- expand.grid(p3, aa)
    pept <- apply(p4, 1, function(x) paste(x, collapse = ""))
    gram$peptide <- pept
    gram$frequency <- vector()
    return(gram)
  }
}


## ----------------------------------------------------------- ##
#              ngraMatrix(data, n, silent)                      #
## ----------------------------------------------------------- ##
#' Compute n-Gram Frequencies Dataframe
#' @description Computes the n-gram frequencies dataframe for the protein and species provides.
#' @usage ngraMatrix(data, n = 2, silent = FALSE)
#' @param data a dataframe with as many columns as species and one row per orthologous protein. The rows and columns must be named accordingly.
#' @param silent logical, set to FALSE to avoid loneliness.
#' @param n a positive integer between 1 and 4.
#' @details The argument prot can be obtained using orth() and orth.seq().
#' @return A list with two dataframes. The first one with nsp * npr columns (nsp: number of species, npr: number of proteins per species) and npe rows (npe: number of peptides, 20 for n = 1, 400 for n = 2, 8000 for n = 3 and 160000 for n = 4). The entries of the dataframe are the number of times that the indicated peptide has been counted in the given protein. Orthologous proteins are in consecutive columns, thus the first nsp columns are the orthologous of protein 1 and so on. The second dataframe contains the Species Vector Sums (each vector describes one species).
#' @references Stuart et al. Bioinformatics 2002; 18:100-108.
#' @seealso ngram(), svdgram()
#' @examples ngraMatrix(bovids[,1:3], n = 2)
#' @export

ngraMatrix <- function(data, n = 2, silent = FALSE){
  aa <- c("A","C","D","E","F","G","H","I","K","L","M",
          "N","P","Q","R","S","T","V","W","Y", "X")
  nsp <- length(data) # number of species
  npr <- nrow(data) # number of proteins
  npe <- 20^n # number of peptides
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
      t <- ngram(data[i,j], n = n)
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
#' a <- ngraMatrix(bovids[, 1:4], n = 2)[[2]][, -1]
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

## ----------------------------------------------------------- ##
#              vtree(matrix, outgroup)                    #
## ----------------------------------------------------------- ##
#' Build a Tree When Wpecies Are Encoded by n-Dim Vectors
#' @description Builds a tree when species are encoded by n-dim vectors.
#' @usage vtree(matrix, outgroup = 'any')
#' @param matrix either a dataframe or matrix where each column represents an OTU.
#' @param outgroup when a rooted tree is desired, it indicates the species to be used as outgroup.
#' @details The method is based on a distance matrix obtained after converting the cos between vector (similarity measurement) in a dissimilarity measurement.
#' @return A list with two objects, the first one is an inter-species distance matrix. The second one is an object of class 'phylo'.
#' @seealso svdgram
#' @examples
#' data(bovids)
#' mymatrix <- ngraMatrix(bovids[, 6:11], n = 2)[[2]][, 2:7]
#' vtree(mymatrix, outgroup = "Pseudoryx_nghetinhensis")
#' @importFrom ape nj
#' @importFrom ape root
#' @importFrom ape plot.phylo
#' @export

vtree <- function(matrix, outgroup = "any"){

  cosDist <- vcos(matrix, silent = TRUE, digits = 6)
  d <- vdis(cosDist)
  d[is.na(d)] <- 0
  d <- d + t(d)
  t <- ape::nj(d)
  if (outgroup[1] != "any"){
    t <- ape::root(t, outgroup = outgroup)
  }
  ape::plot.phylo(t, use.edge.length = FALSE, cex = 0.5)

  return(list(d, t))
}
