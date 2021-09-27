## --------- metrics.R ------------- ##
#                                     #
#      metrics                        #
#      vcos                           #
#      cos2dis                        #
#      vtree                          #
## --------------------------------- ##

## ----------------------------------------------------------- ##
#             metrics(vset, method = 'euclidean', p = 2)        #
## ----------------------------------------------------------- ##
#' Pairwise Vector Dissimilarities
#' @description Computes the dissimilarity between n-dimensional vectors.
#' @usage metrics(vset, method = 'euclidean', p = 2)
#' @param vset matrix (n x m) where each column is a n-dimensional vector.
#' @param method a character string indicating the distance/dissimilarity method to be used (see details).
#' @param p power of the Minkowski distance. This parameter is only relevant if the method 'minkowski' has been selected.
#' @details Although many of the offered methods compute a proper distance, that is not always the case. For instance, for a non null vector, v, the 'cosine' method gives d(v, 2v) = 0, violating the coincidence axiom. For that reason we prefer to use the term dissimilarity instead of distance. The methods offered can be grouped into families.
#'
#' L_p family:
#' ---------------------------------------------------------------
#' ('euclidean', 'manhattan', 'minkowski', 'chebyshev')
#'
#' Euclidean = sqrt( sum | P_i - Q_i |^2)
#'
#' Manhattan = sum | P_i - Q_i |
#'
#' Minkowski = ( sum| P_i - Q_i |^p)^1/p
#'
#' Chebyshev = max | P_i - Q_i |
#'
#' L_1 family:
#' ---------------------------------------------------------------
#' ('sorensen', 'soergel', 'lorentzian', 'kulczynski', 'canberra')
#'
#' Sorensen = sum | P_i - Q_i | / sum (P_i + Q_i)
#'
#' Soergel = sum | P_i - Q_i | / sum max(P_i , Q_i)
#'
#' Lorentzian = sum ln(1 + | P_i - Q_i |)
#'
#' Kulczynski = sum | P_i - Q_i | / sum min(P_i , Q_i)
#'
#' Canberra = sum | P_i - Q_i | / (P_i + Q_i)
#'
#' Intersection family:
#' ---------------------------------------------------------------
#' ('non-intersection', 'wavehedges', 'czekanowski', 'motyka')
#'
#' Non-intersection = 1 - sum min(P_i , Q_i)
#'
#' Wave-Hedges = sum | P_i - Q_i | / max(P_i , Q_i)
#'
#' Czekanowski = sum | P_i - Q_i | / sum | P_i + Q_i |
#'
#' Motyka = sum max(P_i , Q_i) / sum (P_i , Q_i)
#'
#' Inner product family:
#' ---------------------------------------------------------------
#' ('cosine', 'jaccard')
#'
#' Cosine = - ln(0.5 (1 +  (P_i Q_i) / sqrt(sum P_i^2) sqrt(sum Q_i^2)))
#'
#' Jaccard = 1 - sum (P_i Q_i) / (sum P_i^2 + sum Q_i^2 - sum (P_i Q_i))
#'
#' Squared-chord family:
#' ---------------------------------------------------------------
#' ('bhattacharyya', 'squared_chord')
#'
#' Bhattacharyya = - ln sum sqrt(P_i Q_i)
#'
#' Squared-chord = sum ( sqrt(P_i) - sqrt(Q_i) )^2
#'
#' Squared Chi family:
#' ---------------------------------------------------------------
#' ('squared_chi')
#'
#' Squared-Chi = sum ( (P_i - Q_i )^2 / (P_i + Q_i) )
#'
#' Shannon's entropy family:
#' ---------------------------------------------------------------
#' ('kullback-leibler', 'jeffreys', 'jensen-shannon', 'jensen_difference')
#'
#' Kullback-Leibler = sum P_i * log(P_i / Q_i)
#'
#' Jeffreys = sum (P_i - Q_i) * log(P_i / Q_i)
#'
#' Jensen-Shannon = 0.5(sum P_i ln(2P_i / (P_i + Q_i)) + sum Q_i ln(2Q_i / (P_i + Q_i)))
#'
#' Jensen difference = sum (0.5(P_i log(P_i) + Q_i log(Q_i)) - 0.5(P_i + Q_i) ln(0.5(P_i + Q_i))
#'
#' Mismatch family:
#' ---------------------------------------------------------------
#' ('hamming', 'mismatch', 'mismatchZero', 'binary')
#'
#' Hamming = (# coordinates where P_i != Q_i) / n
#'
#' Mismatch = # coordinates where that P_i != Q_i
#'
#' MismatchZero = Same as mismatch but after removing the coordinates where both vectors have zero.
#'
#' Binary = (# coordinates where a vector has 0 and the other has a non-zero value) / n.
#'
#' Combinations family:
#' ---------------------------------------------------------------
#' ('taneja', 'kumar-johnson', 'avg')
#'
#' Taneja = sum ( P_i + Q_i / 2) log( P_i + Q_i / ( 2 sqrt( P_i * Q_i)) )
#'
#' Kumar-Johnson = sum (P_i^2 - Q_i^2)^2 / 2 (P_i Q_i)^1.5
#'
#' Avg = 0.5 (sum | P_i - Q_i| + max{ | P_i - Q_i |})
#'
#' @return A matrix with the computed dissimilarity values.
#' @references Sung-Hyuk Cha (2007). International Journal of Mathematical Models and Methods in Applied Sciences. Issue 4, vol. 1
#' @references Luczac et al. (2019). Briefings in Bioinformatics 20: 1222-1237.
#' @references https://r-snippets.readthedocs.io/en/latest/real_analysis/metrics.html
#' @seealso vcos(), vdis()
#' @examples metrics(matrix(1:9, ncol =3), 'cosine')
#' @importFrom stats dist
#' @importFrom philentropy distance
#' @export

metrics <- function(vset, method = 'euclidean', p = 2){


  if (method %in% c('euclidean', 'chebyshev', 'manhattan',
                    'minkowski', 'canberra', 'binary')){
    if (method == 'chebyshev') {method <- 'maximum'} # Renamed.

    d <- as.matrix(stats::dist(t(vset), method = method, p = p))

  } else if (method %in% c('sorensen', 'soergel', 'kulczynski', 'avg',
                           'lorentzian', 'wavehedges', 'jaccard', 'czekanowski',
                           'squared_chi', 'motyka', 'squared_chord')){
    if (method == 'kulczynski'){ method <- 'kulczynski_d'} # Renamed.

    d <- philentropy::distance(t(vset), method = method, use.row.names = TRUE)

  } else if (method %in% c('non-intersection','bhattacharyya',
                           'hellinger', 'kullback-leibler', 'jeffreys',
                           'jensen-shannon', 'taneja', 'kumar-johnson',
                           'jensen_difference')){
    # counts are converted to frequencies before computation of distances

    d <- philentropy::distance(t(vset), method = method, est.prob = 'empirical',
                               use.row.names = TRUE)

  } else if (method == 'cosine'){
    # the function from EnvNJ package is faster than the function
    # from the philentropy package.

    d <- cos2dis(vcos(vset))
    d[which(is.na(d))] <- 0
    d <- d + t(d)

  } else if (method == 'hamming'){
    # From the Mismatch family

    d <- matrix(NA, nrow = dim(vset)[2], ncol = dim(vset)[2])
    for ( i in 1:((dim(vset)[2])-1) ){
      v <- vset[,i]
      for ( j in (i+1):dim(vset)[2] ){
        w <- vset[,j]
        d[i,j] <- sum(v != w)/dim(vset)[1]
      }
    }
    d[is.na(d)] <- 0
    d <- d + t(d)
    colnames(d) <- rownames(d) <- colnames(vset)

  } else if (method == 'mismatch'){

      d <- matrix(NA, nrow = dim(vset)[2], ncol = dim(vset)[2])
      for ( i in 1:((dim(vset)[2])-1) ){
        v <- vset[,i]
        for ( j in (i+1):dim(vset)[2] ){
          w <- vset[,j]
          d[i,j] <- sum(v != w)
        }
      }
      d[is.na(d)] <- 0
      d <- d + t(d)
      colnames(d) <- rownames(d) <- colnames(vset)

  } else if (method == 'mismatchZero'){
      # Same as mismatch but after removing the
      # coordinates where both vectors have zero.

      d <- matrix(NA, nrow = dim(vset)[2], ncol = dim(vset)[2])
      for ( i in 1:((dim(vset)[2])-1) ){
        v <- vset[,i]
        for ( j in (i+1):dim(vset)[2] ){
          w <- vset[,j]
          Nzeros <-length(which(v == w & v ==0))
          d[i,j] <- sum(v != w)/(dim(vset)[1] - Nzeros)
        }
      }
      d[is.na(d)] <- 0
      d <- d + t(d)
      colnames(d) <- rownames(d) <- colnames(vset)
  }

  return(d)
}


## ----------------------------------------------------------- ##
#               vcos(vectors, silent, digits = 3)               #
## ----------------------------------------------------------- ##
#' Compute Pairwise Cosines of the Angles Between Vectors
#' @description Computes pairwise cosines of the angles between vectors.
#' @usage vcos(vectors, silent = TRUE, digits = 3)
#' @param vectors a named list (or dataframe) containing n-dimensional vectors.
#' @param silent logical, set to FALSE to avoid loneliness.
#' @param digits integer indicating the number of decimal places.
#' @details Cosines are standard measure of vector similarity. If the angle between two vectors in n-dimensional space is small, then the individual elements of their vectors must be very similar to each other in value, and the calculated cosine derived from these values is near one. If the vectors point in opposite directions, then the individual elements of their vectors must be very dissimilar in value, an the calculated cosine is near minus one.
#' @return A triangular matrix with the cosines of the angles formed between the  given vectors.
#' @seealso vdis()
#' @examples vcos(otu.space(bovids[, 1:4]))
#' @export

vcos <- function(vectors, silent = TRUE, digits = 3){

  vectors <- as.data.frame(vectors)

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


## ---------------------------------------------------- ##
##                   cos2dis()                          ##
## ---------------------------------------------------- ##
#' Convert Cosines Between Vectors into Pairwise Dissimilarities
#' @description Converts cosines values into dissimilarities values.
#' @usage cos2dis(cos)
#' @param cos a square upper triangular matrix where cos(i,j) is the cosine between the vector i and j.
#' @details Cosines are standard measure of vector similarity, and can be converted into distance by dij = -log( (1 + cos(i,j) )/2).
#' @return A triangular matrix with the distances.
#' @seealso vcos()
#' @examples
#' data(bovids)
#' vectors = otu.space(bovids[, 7:11])
#' cosData = vcos(vectors)
#' disData = cos2dis(cosData)
#' @export
#'
cos2dis <- function(cos){

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
## ----------------------------------------------------------- ##
#              vtree(matrix, outgroup)                    #
## ----------------------------------------------------------- ##
#' Build a Tree When Species Are Encoded by n-Dim Vectors
#' @description Builds a tree when species are encoded by n-dim vectors.
#' @usage vtree(matrix, outgroup = 'any')
#' @param matrix either a dataframe or matrix where each column represents an OTU.
#' @param outgroup when a rooted tree is desired, it indicates the species to be used as outgroup.
#' @details The method is based on a distance matrix obtained after converting the cos between vector (similarity measurement) in a dissimilarity measurement.
#' @return A list with two objects, the first one is an inter-species distance matrix. The second one is an object of class 'phylo'.
#' @seealso svdgram
#' @examples
#' data(bovids)
#' mymatrix <- ngraMatrix(bovids[, 6:11], k = 2)[[2]][, 2:7]
#' vtree(mymatrix, outgroup = "Pseudoryx_nghetinhensis")
#' @importFrom ape nj
#' @importFrom ape root
#' @importFrom ape plot.phylo
#' @export

vtree <- function(matrix, outgroup = "any"){

  cosDist <- vcos(matrix, silent = TRUE, digits = 6)
  d <- cos2dis(cosDist)
  d[is.na(d)] <- 0
  d <- d + t(d)
  t <- ape::nj(d)
  if (outgroup[1] != "any"){
    t <- ape::root(t, outgroup = outgroup)
  }
  ape::plot.phylo(t, use.edge.length = FALSE, cex = 0.5)

  return(list(d, t))
}
