## ----------- ncdnj.R ------------- ##
#                                     #
#    ncd                              #
#    ncdnj                            #
#                                     #
## --------------------------------- ##

## ----------------------------------------------------------- ##
#                      ncd(seq1, seq2)                          #
## ----------------------------------------------------------- ##
#' Compute Normalized Compression Distances
#' @description Computes normalized compression distances.
#' @usage ncd(seq1, seq2)
#' @param seq1 character string indicating the path to the first fasta file to be analyzed.
#' @param seq2 character string indicating the path to the second fasta file to be analyzed.
#' @details The two fasta files must be in the working directory. This function use zpaq to compress files. Thus, the zpaq software must be installed on your system and in the search path for executables if you wish to use this function.
#' NCD = (Z(xy) - min(Z(x), Z(y))) / max(Z(x), Z(y))
#' Where Z(x), Z(y) and Z(xy) are the lengths of the compressed versions of seq1, seq2 and the concatenated sequences 1 and 2, respectively.
#' @return A non-negative real value reflecting the dissimilarity between seq1 and seq2.
#' @seealso ncdnj()
#' @examples \donttest{try(ncd(seq1 = "./A.fasta", seq2 = "./B.fasta"))}
#' @export

ncd <- function(seq1, seq2){

  t1 <- strsplit(seq1, split = ".", fixed = TRUE)[[1]][1]
  t2 <- strsplit(seq2, split = ".", fixed = TRUE)[[1]][1]
  # ext <- strsplit(seq1, split = ".", fixed = TRUE)[[1]][2]

  ## --- seq1 and seq2 concatenation into s1s2 file
  con <- file(paste(seq1, sep = ""))
  open(con, "r")
  s1 <- readLines(con)
  close(con)

  con <- file(paste(seq2, sep = ""))
  open(con, "r")
  s2 <- readLines(con)
  close(con)

  s1s2 <- paste(s1, s2, sep = "")
  con <- file("s1s2.fasta", "w")
  writeLines(s1s2, con)
  close(con)

  ## --- Compresing seq1, seq2 and s1s2 files
  temp1 <- paste(t1, ".zpaq", sep = "")
  com <- paste("zpaq a ", temp1, " ", seq1, sep = "")
  suppressWarnings(system(com))
  Zseq1 <- file.size(temp1)

  temp2 <- paste(t2, ".zpaq", sep = "")
  com <- paste("zpaq a ", temp2, " ", seq2, sep = "")
  suppressWarnings(system(com))
  Zseq2 <- file.size(temp2)


  suppressWarnings(system("zpaq a s1s2.zpaq s1s2.fasta"))
  Zseq1seq2 <- file.size("s1s2.zpaq")

  ## --- Cleaning up
  if (file.exists('s1s2.zpaq')){
    file.remove('s1s2.zpaq')
  }
  if (file.exists(temp1)){
    file.remove(temp1)
  }
  if (file.exists(temp2)){
    file.remove(temp2)
  }
  if (file.exists("./s1s2.fasta")){
    file.remove("./s1s2.fasta")
  }


  ## --- Computing NCD
  ncd <- (Zseq1seq2 - min(Zseq1, Zseq2))/max(Zseq1, Zseq2)

  return(ncd)
}

## ----------------------------------------------------------- ##
#                            ncdnj(wd)                          #
## ----------------------------------------------------------- ##
#' Compute a Distance Matrix Using Normalized Compression Distance
#' @description Computes a distance matrix using normalized compression distance.
#' @usage ncdnj(wd)
#' @param wd character string indicating the path to the directory where the input files can be found (see details).
#' @details The input files, which must be found at the wd provided, consist of a file named 'list.txt' containing the names of the fasta files to be analyzed (one per line). The referred fasta files also must be found at the provided wd. This function use zpaq to compress files. Thus, the zpaq software must be installed on your system and in the search path for executables if you wish to use this function.
#' @return A list where the first element is a symmetric distance matrix and the second one is a phylogenetic tree build using NJ.
#' @seealso ncd()
#' @examples \donttest{try(ncdnj("./data_t"))}
#' @importFrom ape nj
#' @export

ncdnj <- function(wd){

  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd(wd)

  l <- unname(unlist(read.csv("./list.txt", header = FALSE)))
  sp <- unlist(lapply(l, function(x) gsub(".fasta", "", x)))

  d <- matrix(rep(NA, length(l)^2), ncol = length(l))
  colnames(d) <- rownames(d) <- sp

  for (i in 1:(length(l) - 1)){
    for (j in (i+1):length(l)){
      d[i,j] <- ncd(seq1 = l[i], seq2 = l[j])
    }
  }
  d[which(is.na(d))] <- 0
  d <- d + t(d)
  if (dim(d)[1] > 3){
    tree <- ape::nj(d)
  } else {
    tree <- NULL
  }

  return(list(d, tree))
}

