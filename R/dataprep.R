## ---------- dataprep.R ----------- ##
#                                     #
#    fastaconc                        #
#    df2fasta                         #
#    d.phy2df                         #
#                                     #
## --------------------------------- ##

## ----------------------------------------------------------- ##
#           fastaconc(otus, inputdir,out.file)                  #
## ----------------------------------------------------------- ##
#' Concatenate Fasta Files in a Single Multispecies Fasta File
#' @description Concatenate fasta files from different species in a single multispecies fasta file.
#' @usage fastaconc(otus, inputdir = ".", out.file = "./concatenated_multispecies.fasta")
#' @param otus a character vector giving the otus' names.
#' @param inputdir path to the directory containing the individual fasta files.
#' @param out.file path and name of output file.
#' @details When we have fasta files (extension should be '.fasta'), each one for a species containing different sequences of the given species, this function concatenate the different sequences of the same species and writes it as a single sequence in a single multispecies fasta file. If the individual fasta files are found in the working directory, the inputdir argument don't need to be passed. The names of the individual fasta files must match the otus' names.
#' @return A single multispecies fasta file with the sequences of each species spliced in a single sequence.
#' @seealso df2fasta()
#' @examples \dontrun{fastaconc(otus = c('Glis_glis', 'Ovis_aries', 'Sus_scrofa'))}
#' @importFrom seqinr read.fasta
#' @importFrom seqinr write.fasta
#' @export

fastaconc <- function(otus, inputdir = ".", out.file = "./concatenated_multispecies.fasta"){

  seqs <- character(length(otus))
  for (i in 1:length(otus)){
    f <- unlist(seqinr::read.fasta(paste(inputdir, "/", otus[i], ".fasta", sep = "")))
    seqs[i] <- paste(toupper(f), collapse = "")
  }

  seqinr::write.fasta(sequences = as.list(seqs),
                      names = otus,
                      file.out = out.file,
                      as.string = TRUE)

  print(paste("Work finished. Fasta file saved at", out.file))
}


## ---------------------------------------------------- ##
##         df2fasta(df, out.file)                       ##
## ---------------------------------------------------- ##
#' Convert Dataframe into Fasta File
#' @description Converts a dataframe into a fasta file.
#' @usage df2fasta(df, out.file)
#' @param df a named (both rows and cols) dataframe (see details).
#' @param out.file path and name of output file.
#' @details The format of the df should be as follows. Each row represents a protein sequence and each column a species.
#' @return A fasta file that is saved in the specified path.
#' @seealso fastaconc()
#' @examples \dontrun{df2fasta(df = bovids, out.file = "./example.fasta")}
#' @importFrom seqinr write.fasta
#' @importFrom seqinr s2c
#' @export

df2fasta <- function(df, out.file){

  # protname <- rownames(df)
  otus <- names(df)

  seqs <- unlist(lapply(1:length(otus), function(i) paste(df[,i], collapse = "")))
  # seqs <- lapply(seqs, seqinr::s2c)

  seqinr::write.fasta(sequences = as.list(seqs),
                        names = otus,
                        file.out = out.file,
                        as.string = TRUE,)


  print(paste("Work finished. Files saved at", out.file))
}


## ---------------------------------------------------- ##
##          d.phy2df(phyfile, as)                       ##
## ---------------------------------------------------- ##
#' Convert a Phylip Distance Matrix into a DataFrame
#' @description Converts a phylip distance matrix into a either an R dataFrame or matrix.
#' @usage d.phy2df(phyfile, as = 'matrix')
#' @param phyfile path to the file containing the distances in phylip format.
#' @param as class of the output R data. It should be either 'dataframe' or 'matrix'.
#' @return Either a dataframe or a matrix containing the distances read from the phy file.
#' @seealso d.df2pny()
#' @examples \dontrun{d.phy2df(phyfile = "./data_t/d_dummy.txt")}
#' @importFrom utils read.csv
#' @export
d.phy2df <- function(phyfile, as = 'matrix'){
  d <- read.csv(phyfile, sep = "", header = FALSE, skip = 1)
  names <- d[,1]
  d <- d[,-1]
  names(d) <- names
  rownames(d) <- names
  if (as == 'matrix'){
    d <- as.matrix(d)
  }
  return(d)
}

