## ------- envcoding.R ------------ ##
#                                    #
#   env.extract                      #
#   env.matrices                     #
#   env.sp                           #
#   otu.vector                       #
#   otu.space                        #
#                                    #
## -------------------------------- ##


## ------------------------------------------------------------------------------- ##
#                   env.extract <- function(seq, c, r)                              #
## ------------------------------------------------------------------------------- ##
#' Sequence Environment Around a Given Position
#' @description Extracts the sequence environment around a given position.
#' @usage env.extract(seq, c, r)
#' @param seq a string protein sequence.
#' @param c center of the environment.
#' @param r radius of the environment.
#' @return Returns a  a strings with the extracted environment sequence.
#' @examples env.extract('ARGGQQQCATSYV', c = 8,  r = 2)
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @seealso env.matrices(), env.sp()
#' @importFrom bio3d get.seq
#' @export

env.extract <- function(seq,  c, r){

  ## ---- Checking the request's coherence ---- ##
  if (c > nchar(seq)){
    stop("The requested center is not found in the provided sequence")
  }

  ## -- Ancillary function to isolate environment -- ##
  extract <- function(seq, c, r){
    envL <- substring(seq, c-r, c-1)
    if (nchar(envL) < r){
      X <- paste(rep('X', r-nchar(envL)), collapse = "")
      envL <- paste(X, envL, sep = "")
    }

    envR <- substring(seq, c+1, c+r)
    if (nchar(envR) < r){
      X <- paste(rep('X', r-nchar(envR)), collapse = "")
      envR <- paste(envR, X, sep = "")
    }

    envC <- tolower(substring(seq, c, c))
    env <- paste(envL, envC, envR, sep = "")
    return(env)
  }

  return(extract(seq, c, r))
}


## ---------------------------------------------------------------- ##
#               env.matrices <- function(env)                        #
## ---------------------------------------------------------------- ##
#' Environment Matrices
#' @description Provides the frequencies of each amino acid within the environment.
#' @usage env.matrices(env)
#' @param env a character string vector containing the environments.
#' @return Returns a list of two dataframes. The first, shown the environment in matrix form. The second provides the frequencies of each amino acid within the environments.
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @examples env.matrices(c('ANQRmCTPQ', 'LYPPmQTPC', 'XXGSmSGXX'))
#' @seealso env.extract(), otu.vector()
#' @export

env.matrices <- function(env){

  ## --------- Checking input data --------- ##
  env <- env[!is.na(env)] # remove NA if needed
  L <- unique(nchar(env)) # environment's length
  if (length(L) != 1){
    stop("The length of sequence environments is not homogeneous in your input data")
  }
  # # Remove env with non-canonical amino acids:
  env <- gsub("[^A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
               a, c, d, e, f, g, h, i, k, l, m, n, p, q, r, s, t, v, w, y]",
              "", env)
  if (length(which(nchar(env) != L)) > 0){
    env <- env[- which(nchar(env) != L)]
  }

  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X")
  n_env <- length(env) # number of environments being analyzed

  r <- L %/% 2 # environment's radius

  ## ---- Dealing with the amino acids DF ---- ##
  g <- mapply(strsplit,  env, "")
  aaDF <- data.frame(do.call(rbind, g))
  names(aaDF) <- -r:r

  if (nrow(aaDF) == 1){ # working copy when a single environment seq
    aaDF_ <- as.data.frame(t(apply(aaDF, 2, toupper)))
    rownames(aaDF_) <- rownames(aaDF)
  } else { # working copy when several environment seqs
    aaDF_ <- as.data.frame(apply(aaDF, 2, toupper))
  }
  aaDF_ <- as.data.frame(lapply(aaDF_, factor, levels = aa))

  ## ---- Dealing with the frequencies DF ---- ##
  fDF <- as.data.frame(matrix(rep(NA, L*21), ncol = L))
  names(fDF) <- -r:r
  rownames(fDF) <- aa

  for (i in 1:L){
    fDF[,i] <- table(aaDF_[,i])
  }
  t <- toupper(aaDF[,r+1][1]) # central residue
  if (fDF[which(aa == t), r+1] != n_env){
    warn <- paste("The amino acid found at the central position of the environment \n",
                  "may not be always the same in your input data", sep = "")
    warning(warn)
  }
  attr(fDF, 'number_sequences_analyzed') <- n_env
  output <- list(aaDF, fDF)

  return(output)
}


## ----------------------------------------------------------- ##
#               env.sp(data, sp, r, aa, silent)                 #
## ----------------------------------------------------------- ##
#' Extract the Sequence Environments
#' @description Extracts the sequence environments around the selected amino acid(s) in the chosen species.
#' @usage env.sp(data, sp, r = 10, aa = 'all', silent = TRUE)
#' @param data input data must be a dataframe (see details).
#' @param sp the species of interest (it should be named as in the input dataframe).
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) which environments are going to be extracted.
#' @param silent logical. When FALSE the program progress is reported to alleviate loneliness.
#' @details Input data must be a dataframe where each row corresponds to an individual protein. The columns contain the sequence of the protein corresponding to the row in each species. Therefore, the columns' names of this dataframe must be coherent with the names of the OTUs being analyzed.
#' @return A list of environment sequences. Each element from the list is a vector with the environment sequences around an amino acid. So, the length list is the same as the length of aa.
#' @seealso otu.vector(), otu.space()
#' @examples
#' data(bovids)
#' env.sp(data = bovids, sp = "Bos_taurus", r = 2)
#' @importFrom stringr str_which
#' @export

env.sp <- function(data, sp, r = 10, aa = 'all', silent = TRUE){

  co <- which(names(data) == sp) # column from data corresponding to the relevant species
  separador <- paste(rep("X", r), collapse = "")
  p <- paste(data[, co], collapse = separador) # single sequence

  aminoacidos <-  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  if (aa[1] == 'all'){
    aa <- aminoacidos
  }

  env <- vector(mode = "list", length = length(aa)) # will contain the environment for each aa
  names(env) <- aa
  counter <- 0
  for (a in aa){
    if (!silent){
      print(paste("  -----  ", a, "  -----  ", sp, "  -----"), sep = "")
    }
    counter = counter + 1
    t <- gregexpr(a, p)[[1]]
    if (t[1] != -1){ # Just in case the protein doesn't have the amino acid bein analyzed
      temp <- unlist(lapply(t, function(x) env.extract(p, c = x, r = r)))
      tbr <- stringr::str_which(temp, "X")
      if (length(tbr != 0)){ # There must be r residues at each end of the central aa
        env[[counter]] <- temp[-tbr]
      } else {
        env[[counter]] <- temp
      }
    }
  }
  return(env)
}


## ----------------------------------------------------------- ##
#                 otu.vector(envl, sp, aa, silent)              #
## ----------------------------------------------------------- ##
#' Convert a Set of Sequence Environments into a Vector
#' @description Converts a set of sequence environments into a vector.
#' @usage otu.vector(envl, sp = "", aa = "all", silent = TRUE)
#' @param envl a list containing the sequence environment of a species (as the one returned by the function env.sp()).
#' @param sp character string indicating the species being analyzed.
#' @param aa the amino acid(s) to be used to encoded the species.
#' @param silent logical. When FALSE the program progress is reported to alleviate loneliness.
#' @details The dimension of the vector representing the species will depend on the settings. For instance, if we choose a single amino acid and a radius of 10 for the sequence environment, then we will get a vector of dimension 400 (20 amino acids x 20 positions). If we opt for the 20 amino acids and r = 10, then the vector will be of dimension 8000 (400 for each amino acid * 20 amino acids). Please, note that r is selected in the function env.sp() that will provide the input dataframe for the current function.
#' @return A matrix representing the species. This matrix can be converted into a vector representing the target species just typing as.vector(matrix). Each coordinate is the frequency of a given amino acid at a certain position from the environment (see details).
#' @seealso env.sp(), otu.space()
#' @examples
#' data(bovids)
#' cow = env.sp(bovids, "Bos_taurus")
#' otu.vector(cow)
#' @export

otu.vector <- function(envl, sp = '', aa = "all", silent = TRUE){

  aminoacidos <-  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  ## ----- check that aa are licit
  if (aa[1] != "all"){
    if (sum(!aa %in% aminoacidos) > 0){
      stop(paste("No allowed: ", aa[!aa %in% aminoacidos], "\n"))
    } else {
      aminoacidos <- aa
    }
  }

  pep <- envl[[1]][1]
  npos <- nchar(pep) - 1 # number of position around target amino acid
  center <- (npos/2) + 1 # center (position of target amino acid)
  Vmatrix <- matrix(rep(NA, npos*20 * length(aminoacidos)),
                    ncol = length(aminoacidos))
  colnames(Vmatrix) <- aminoacidos

  for (i in 1:length(aminoacidos)){ # for each amino acid
    aat <- aminoacidos[i]
    if (!silent){
      print(paste("-------- ", aat, " --------"))
    }

    env <- envl[[i]]
    if (length(env) == 0){
      Vmatrix[, i] <- as.vector(rep(0, nrow(Vmatrix)))
      next # When the amino acid is not present in the sequence
    } else {
      A <- env.matrices(env)[[2]]
      ## -- Quality test
      if (length(unique(apply(A, 2,sum))) != 1){
        stop(paste("All columns must add up to the same", sp))
      } else if (unique(apply(A, 2,sum)) != A[which(rownames(A) == aat), center]){
        stop(paste("All columns must add up to number of aa", sp))
      } else if (sum(A[21, ]) != 0){
        stop(paste("X is not an expected character at", sp))
      }

      Abs <- as.matrix(A[1:20, -center]) # remove the last row (X character) and the central col
      Vmatrix[, i] <- as.vector(Abs)
    }
  }

  attr(Vmatrix, "species") <- sp
  attr(Vmatrix, "aa") <- aminoacidos

  return(Vmatrix)
}


## ----------------------------------------------------------- ##
#                 otu.space(data, r, aa, silent)                #
## ----------------------------------------------------------- ##
#' Compute the Matrix Representing the Species Vector Subspace
#' @description Computes the matrix representing the species vector subspace.
#' @usage otu.space(data, r = 10, aa = "all", silent = TRUE)
#' @param data input data must be a dataframe (see details).
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) to be used to encoded the species.
#' @param silent logical. If FALSE, the running progress is reported.
#' @details Input data must be a dataframe where each row corresponds to an individual protein, and each column identifies a species. Therefore, the columns' names of this dataframe must be coherent with the names of the OTUs being analyzed. The dimension of the vector representing each species will depend on the settings. For instance, if we choose a single amino acid and a radius of 10 for the sequence environment, then we will get a vector of dimension 400 (20 amino acids x 20 positions). If we opt for the 20 amino acids and r = 10, then the vector will be of dimension 8000 (400 for each amino acid * 20 amino acids). Please, note that r is selected in the function env.sp() that will provide the input dataframe for the current function.
#' @return A matrix representing the species vector subspace.
#' @seealso env.sp(), otu.vector()
#' @examples
#' data(bovids)
#' otu.space(bovids[, 1:5], r = 2)
#' @export

otu.space <- function(data, r = 10, aa = "all", silent = TRUE){

  species <- names(data)

  mylist <- lapply(species, function(x) env.sp(data = data, x, r = r, silent = silent))
  myspace <- lapply(mylist, function(x) as.vector(otu.vector(envl = x, aa = aa)))
  space <- matrix(unlist(myspace), ncol = length(species), byrow = FALSE)
  colnames(space) <- species

  return(space)
}

