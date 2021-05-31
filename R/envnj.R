## ----------- envnj.R ------------- ##
#                                     #
#    envnj                            #
#    vcos                             #
#    vdis                             #
#                                     #
## --------------------------------- ##


## ----------------------------------------------------------- ##
#               envnj(data, r, aa, outgroup)                    #
## ----------------------------------------------------------- ##
#' Build Trees Based on the Environment Around the Indicated Amino Acid(s)
#' @description Builds trees based on the environment around the indicated amino acid(s).
#' @usage envnj(data, r = 10, aa = 'all', outgroup = 'any')
#' @param data input data must be a dataframe where each row corresponds to a protein sequence and each column to a species.
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) to be used to encoded the species.
#' @param outgroup when a rooted tree is desired, it indicates the species to be used as outgroup.
#' @details This function builds alignment-independent phylogenetic trees.
#' @return A list with two objects, the first one is an inter-species distance matrix. The second one is an object of class 'phylo'.
#' @seealso otu.space()
#' @examples \donttest{
#' data(bovids)
#' envnj(bovids[, 7:11], aa = "all", outgroup = "Pseudoryx_nghetinhensis")
#' }
#' @importFrom ape nj
#' @importFrom ape root
#' @importFrom ape plot.phylo
#' @export

envnj <- function(data, r = 10, aa = "all", outgroup = 'any'){
  space <- otu.space(data = data, r = r, aa = aa)
  cosDist <- vcos(space, silent = TRUE, digits = 6)
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

## ----------------------------------------------------------- ##
#               vcos(vectors, silent, digits = 3)               #
## ----------------------------------------------------------- ##
#' Compute Pairwise Cosines of the Angles Between Vectors
#' @description Computes pairwise cosines of the angles between vectors.
#' @usage vcos(vectors, silent = FALSE, digits = 3)
#' @param vectors a named list (or dataframe) containing n-dimensional vectors.
#' @param silent logical, set to FALSE to avoid loneliness.
#' @param digits integer indicating the number of decimal places.
#' @details Cosines are standard measure of vector similarity. If the angle between two vectors in n-dimensional space is small, then the individual elements of their vectors must be very similar to each other in value, and the calculated cosine derived from these values is near one. If the vectors point in opposite directions, then the individual elements of their vectors must be very dissimilar in value, an the calculated cosine is near minus one.
#' @return A triangular matrix with the cosines of the angles formed between the  given vectors.
#' @seealso vdis()
#' @examples vcos(otu.space(bovids[, 1:4]))
#' @export

vcos <- function(vectors, silent = FALSE, digits = 3){

  # Check all the vectors have the same dimension
  n <- unlist(lapply(vectors, length))
  if (length(unique(n)) != 1){
    stop("All the vectors must have the same dimension")
  }

  # Matrix to hold the cosines
  n <- length(n) # number of vectors
  cosines <- matrix(rep(NA, n*n), ncol = n)

  for (i in 1:(n-1)){
    if (!silent){print(i)}
    for (j in (i+1):n){
      t <- crossprod(vectors[[i]], vectors[[j]])
      t <- t / sqrt(crossprod(vectors[[i]]) * crossprod(vectors[[j]]))
      cosines[i,j] <- round(t, digits)
    }
  }
  cosines[is.nan(cosines)] <- 0 # The 0 vector is orthogonal to any other vector
  colnames(cosines) <- names(vectors)
  return(cosines)
}

## ----------------------------------------------------------- ##
#                         vdis(cos)                             #
## ----------------------------------------------------------- ##
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
#' disData = vdis(cosData)
#' @export

vdis <- function(cos){

  # Check that cos is a squere matrix
  if (!is.matrix(cos)){
    stop("The argument must be a matrix")
  } else if (length(unique(dim(cos))) != 1){
    stop("The argument matrix must be square")
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
