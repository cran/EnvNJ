## ------- encoding.R ------------- ##
#                                    #
#   aa.at                            #
#   env.extract                      #
#   env.matrices                     #
#   env.sp                           #
#   otu.vector                       #
#   otu.space                        #
#   aaf                              #
#                                    #
## -------------------------------- ##

## --------------------------------------------------------------- ##
#             aa.at(at, target, uniprot = TRUE)                     #
## --------------------------------------------------------------- ##
#' Residue Found at the Requested Position
#' @description Returns the residue found at the requested position.
#' @usage aa.at(at, target, uniprot = TRUE)
#' @param at the position in the primary structure of the protein.
#' @param target a character string specifying the UniProt ID of the protein of interest or, alternatively, the sequence of that protein.
#' @param uniprot logical, if TRUE the argument 'target' should be an ID.
#' @details Please, note that when uniprot is set to FALSE, target can be the string returned by a suitable function, such as get.seq or other.
#' @return Returns a single character representing the residue found at the indicated position in the indicated protein.
#' @examples \donttest{aa.at(28, 'P01009')}
#' @seealso aa.comp()
#' @importFrom bio3d read.fasta
#' @importFrom bio3d get.seq
#' @export

aa.at <- function(at, target, uniprot = TRUE){

  if (uniprot == TRUE){
    target <- tryCatch(
      {
        suppressWarnings(bio3d::get.seq(target)$ali)
      },
      error = function(cond){
        return(NULL)
      }
    )

    if (is.null(target)){
      message("Sorry, get.seq failed")
      return(NULL)
    }

    if (file.exists("seqs.fasta")){
      file.remove("seqs.fasta")
    }

  } else {
    if (length(target) == 1){
      target <- strsplit(target, split="")[[1]]
    }
  }

  if (at %in% 1:length(target)){
    return(target[at])
  } else {
    message(paste(at , " isn't a valid position for this protein", sep=""))
    return(NULL)
  }
}

## ------------------------------------------------------------------------------- ##
#   env.extract <- function(prot, db = 'none', c, r, ctr = 'none', exclude = c())   #
## ------------------------------------------------------------------------------- ##
#' Sequence Environment Around a Given Position
#' @description Extracts the sequence environment around a given position.
#' @usage env.extract(prot, db = 'none', c, r, ctr = 'none', exclude = c())
#' @param prot either a uniprot id or a string sequence.
#' @param db a character string specifying either 'uniprot' or 'none'.
#' @param c center of the environment.
#' @param r radius of the environment.
#' @param ctr the type of control environment; it must be one of 'random', 'closest', or 'none'.
#' @param exclude a vector containing the positions to be excluded as control.
#' @details The random control returns an environment center at a random position containing the same type or amino acid than the positive environment. The closest control searches for the closest position where such a type of amino acid is found and returns its environment.
#' @return Returns a  list of two strings (environments).
#' @examples env.extract('P01009', db = 'uniprot', 271, 10, ctr = 'random')
#' @references Aledo et al. Sci Rep. 2015; 5: 16955. (PMID: 26597773)
#' @seealso env.matrices(), env.sp()
#' @importFrom bio3d get.seq
#' @export

env.extract <- function(prot, db = 'none', c, r, ctr = 'none', exclude = c()){

  ## ------------ Preparing the sequence ------------ ##
  if (db == 'uniprot'){
    seq <- tryCatch(
      {
        bio3d::get.seq(prot)
      },
      error = function(cond){
        return(NULL)
      }
    )
    if (is.null(seq)){
      message("Sorry, get.seq failed")
      return(NULL)
    }

    seq <- paste(seq$ali, collapse = "")

    if (file.exists("seqs.fasta")){
      file.remove("seqs.fasta")
    }

  } else {
    seq <- prot
  }

  ## ---- Checking for coherence of the request ---- ##
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

  ## ------- Building the output file ------- ##

  ## -- The positive:
  env_pos <- extract(seq, c, r)

  ## -- The control:
  res <- aa.at(c, seq, FALSE)
  all_res <- gregexpr(res, seq)[[1]] # Positions containing the same amino acid finds at 'c'
  other_res <- all_res[-which(all_res == c)] # Positions containing the same amino acid except the positive
  other_res <- setdiff(other_res, exclude) # Exclude some positions if requested
  n <- length(other_res)

  if (n == 0){ # There are no other residues in the protein to be used as control
    env_ctr <- ""
  } else {
    if (ctr == 'random'){
      c_random <- other_res[sample(1:n, 1)] # Random position where there is a residue equal to that find at 'c'
      env_ctr <- extract(seq, c_random, r)
    } else if (ctr == 'closest'){
      closest <- which(abs(other_res - c) == min(abs(other_res - c)))[1] # there may be two equidistant residues
      c_closest <- other_res[closest]
      env_ctr <- extract(seq, c_closest, r)
    } else {
      env_ctr <- ""
    }
  }

  output <- list(env_pos, env_ctr)
  names(output) <- c("Positive", "Control")

  return(output)
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
  # non-canonical amino acid to X:
  env <- gsub("[^A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
               a, c, d, e, f, g, h, i, k, l, m, n, p, q, r, s, t, v, w, y]",
              "X", env)

  aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
          "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "X")
  n_env <- length(env) # number of environments being analyzed

  r <- L %/% 2 # environment's radius

  ## ---- Dealing with the amino acids DF ---- ##
  g <- mapply(strsplit,  env, "")
  aaDF <- data.frame(do.call(rbind, g))
  names(aaDF) <- -r:r
  aaDF_ <- as.data.frame(apply(aaDF, 2, toupper)) # working copy
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
#           env.sp(data, sp, r, aa, remove.init, silent)        #
## ----------------------------------------------------------- ##
#' Extract the Sequence Environments
#' @description Extracts the sequence environments around the selected amino acid(s) in the chosen species.
#' @usage env.sp(data, sp, r = 10, aa = 'all', remove.init = TRUE, silent = TRUE)
#' @param data input data must be a dataframe (see details).
#' @param sp the species of interest (it should be named as in the input dataframe).
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) which environments are going to be extracted.
#' @param remove.init logical. If TRUE, the initiation methionine in each protein is disregarded.
#' @param silent logical. When FALSE the program progress is reported to alleviate loneliness.
#' @details Input data must be a dataframe where each row corresponds to an individual protein. The columns contain the sequence of the protein corresponding to the row in each species. Therefore, the columns' names of this dataframe must be coherent with the names of the OTUs being analyzed.
#' @return A dataframe with as many rows as sequence environment have been extracted. The position of the central residue in the protein is also provided.
#' @seealso otu.vector(), otu.space()
#' @examples
#' data(bovids)
#' env.sp(data = bovids, sp = "Bos_taurus", r = 2)
#' @export

env.sp <- function(data, sp, r = 10, aa = 'all', remove.init = TRUE, silent = TRUE){
  co <- which(names(data) == sp) # column from data corresponding to the relevant species

  p <- data[, co] # set of proteins for selected species
  if (remove.init){ # remove initiation methionine if required
    p <- unlist(lapply(p, function(x) substr(x, start = 2, stop = nchar(x))))
  }

  N <- sum(unlist(lapply(p, nchar))) # number of residues
  df <- as.data.frame(matrix(rep(NA, N * 6), ncol = 6))
  names(df) <- c("n", "prot","pos", "env", "species", "aa")

  n <- 1
  aminoacidos <-  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

  for (a in aminoacidos){
    if (!silent){
      print(paste("------ ", a, " -------"))
    }
    met <- lapply(p, function(x) gregexpr(a, x)[[1]])

    lines2beAdded <- length(which(unlist(met) == -1))
    if (lines2beAdded > 0){
      df[nrow(df):(nrow(df) + lines2beAdded), ] <- NA
    }

    for (i in 1:length(met)){ # for each orthologous protein
      if (!silent){
        print(paste("-- Prot ", a, " :", i, sep = ""))
      }
      l <- length(met[[i]])
      df$prot[n:(n+l-1)] <- i
      df$pos[n:(n+l-1)] <- met[[i]]

      # It may happen that the chosen residue is not present in some proteins:
      tryCatch(
        {
          df$env[n:(n+l-1)] <- suppressMessages(unlist(lapply(met[[i]],
                                                  function(x) env.extract(p[i],
                                                          db = "none", c = x,
                                                          r = r, ctr = "none")$Positive)))
        },
        error = function(cond){
          return(NA)
        },
        warning = function(w) conditionMessage(w)
      )

      df$species[n:(n+l-1)] <- sp
      df$aa[n:(n+l-1)] <- a
      n <- n + l
    }
  }

  ## ----- Clean df
  df <- df[df$pos != -1, ]
  emptyrows <- which(is.na(df$pos)) ## ------------------------- Added after gorilla
  if (length(emptyrows) > 0){       ## -------------------------
    df <- df[-emptyrows, ]          ## -------------------------
  }                                 ## -------------------------

  # There must be r residues at each end of the central aa:
  tobedisregarded <- c()
  for (i in 1:nrow(df)){
    if (gregexpr("X", df$env[i])[[1]][1] != -1){
      tobedisregarded <- c(tobedisregarded, i)
    }
  }
  df <- df[-tobedisregarded, ]
  df$n <- 1:nrow(df)

  if (aa[1] != 'all'){
    df <- df[which(df$aa %in% aa), ]
  }

  return(df)
}

## ----------------------------------------------------------- ##
#                      otu.vector(df, aa)                       #
## ----------------------------------------------------------- ##
#' Convert a Set of Sequence Environments into a Vector
#' @description Converts a set of sequence environments into a vector.
#' @usage otu.vector(df, aa = "all")
#' @param df a dataframe containing the sequence environment of a species (as the one returned by the function env.sp()).
#' @param aa the amino acid(s) to be used to encoded the species.
#' @details The dimension of the vector representing the species will depend on the settings. For instance, if we choose a single amino acid and a radius of 10 for the sequence environment, then we will get a vector of dimension 400 (20 amino acids x 20 positions). If we opt for the 20 amino acids and r = 10, then the vector will be of dimension 8000 (400 for each amino acid * 20 amino acids). Please, note that r is selected in the function env.sp() that will provide the input dataframe for the current function.
#' @return A matrix representing the species. This matrix can be converted into a vector representing the target species just typing as.vector(matrix). Each coordinate is the frequency of a given amino acid at a certain position from the environment (see details).
#' @seealso env.sp(), otu.space()
#' @examples
#' data(bovids)
#' cow = env.sp(bovids, "Bos_taurus")
#' otu.vector(cow)
#' @export

otu.vector <- function(df, aa = "all"){

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

  df <- df[which(df$aa %in% aminoacidos), ]
  pep <- df$env[1]
  npos <- nchar(pep) - 1 # number of position around target amino acid
  center <- (npos/2) + 1 # center (position of target amino acid)
  Vmatrix <- matrix(rep(NA, npos*20 * length(aminoacidos)),
                    ncol = length(aminoacidos))
  colnames(Vmatrix) <- aminoacidos

  for (i in 1:length(aminoacidos)){ # for each amino acid
    aat <- aminoacidos[i]
    print(paste("-------- ", aat, " --------"))
    env <- df$env[which(df$aa == aat)]
    A <- env.matrices(env)[[2]]
    ## -- Quality test
    if (length(unique(apply(A, 2,sum))) != 1){
      stop(paste("All columns must add up to the same", df$sp[1]))
    } else if (unique(apply(A, 2,sum)) != A[which(rownames(A) == aat), center]){
      stop(paste("All columns must add up to number of aa", df$sp[1]))
    } else if (sum(A[21, ]) != 0){
      stop(paste("X is not an expected character at ", df$sp[1]))
    }

    Abs <- as.matrix(A[1:20, -center]) # remove the last row (X character) and the central col
    Vmatrix[, i] <- as.vector(Abs)
  }

  attr(Vmatrix, "species") <- df$species[1]
  attr(Vmatrix, "aa") <- aminoacidos

  return(Vmatrix)
}

## ----------------------------------------------------------- ##
#              otu.space(data, r, aa, remove.init)              #
## ----------------------------------------------------------- ##
#' Compute the Matrix Representing the Species Vector Subspace
#' @description Computes the matrix representing the species vector subspace.
#' @usage otu.space(data, r = 10, aa = "all", remove.init = TRUE)
#' @param data input data must be a dataframe (see details).
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param aa the amino acid(s) to be used to encoded the species.
#' @param remove.init logical. If TRUE, the initiation methionine in each protein is disregarded.
#' @details Input data must be a dataframe where each row corresponds to an individual protein, and each column identifies a species. Therefore, the columns' names of this dataframe must be coherent with the names of the OTUs being analyzed. The dimension of the vector representing each species will depend on the settings. For instance, if we choose a single amino acid and a radius of 10 for the sequence environment, then we will get a vector of dimension 400 (20 amino acids x 20 positions). If we opt for the 20 amino acids and r = 10, then the vector will be of dimension 8000 (400 for each amino acid * 20 amino acids). Please, note that r is selected in the function env.sp() that will provide the input dataframe for the current function.
#' @return A matrix representing the species vector subspace.
#' @seealso env.sp(), otu.vector()
#' @examples
#' data(bovids)
#' otu.space(bovids[, 1:5], r = 2)
#' @export

otu.space <- function(data, r = 10, aa = "all", remove.init = TRUE){

  species <- names(data)
  space <- as.data.frame(matrix(rep(NA, 8000*length(species)), ncol = length(species)))
  names(space) <- species

  for (i in 1:length(species)){
    sp <- species[i]
    print(paste("-------- ", sp, " ----------"))
    df <- env.sp(data = data, sp = sp, r = r, remove.init = remove.init )
    if (i == 1){
      v <- as.vector(otu.vector(df = df, aa = aa))
      space <- space[1:length(v), ]
      space[, i] <- v
    }
    space[, i] <- as.vector(otu.vector(df = df, aa = aa))
  }

  return(space)
}

## ----------------------------------------------------------- ##
#                          aaf(data)                            #
## ----------------------------------------------------------- ##
#' Compute the Frequency of Each Amino Acid in Each Species
#' @description Computes the frequency of each amino acid in each species.
#' @usage aaf(data)
#' @param data input data must be a dataframe (see details).
#' @details Input data must be a dataframe where each row corresponds to an individual protein, and each column identifies a species. Therefore, the columns' names of this dataframe must be coherent with the names of the OTUs being analyzed.
#' @return A dataframe providing amino acid frequencies en the set of species. Rows correspond amino acids and columns to species.
#' @seealso env.sp(), otu.vector(), otu.space()
#' @examples
#' data(bovids)
#' aaf(bovids)
#' @importFrom graphics image
#' @export

aaf <- function(data){

  aminoacidos <-  c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  species <- names(data)
  output <- data.frame(matrix(rep(NA, 20 * length(species)), ncol = length(species)))
  names(output) <- species
  rownames(output) <- aminoacidos

  splicedprot <- apply(data, 2, function(x) paste(x, collapse = ""))
  for (i in 1:length(species)){
    output[, i] <- aa.comp(splicedprot[i], uniprot = FALSE)$frequency
  }

  rotate <- function(x) t(apply(x, 2, rev))
  image(rotate(output))

  return(output)
}
