## ----------- envfas.R ------------ ##
#                                     #
#    envfascpp                        #
#    vect2tree                        #
#                                     #
## --------------------------------- ##

## ---------------------------------------------------- ##
##                    envfascpp()                       ##
## ---------------------------------------------------- ##
#' Convert Fasta Files into Environment Vectors
#' @description Converts fasta files into environment vectors
#' @usage envfascpp(path, r = 10, exefile, outfile)
#' @param path path to the folder that contain a the file 'list.txt', which contains the names of all the fasta files to be analyzed (one per line). The fasta files must be in the same path.
#' @param r a positive integer indicating the radius of the sequence segment considered as environment.
#' @param exefile path to the directory containing the envfas executable.
#' @param outfile path to, and name of, output file.
#' @details This function builds and write vectors representing the species' proteome. To use this function, you must first download the source code envfas.cpp (https://bitbucket.org/jcaledo/envnj/src/master/AncillaryCode/envfas.cpp) and compile it.
#' @return A txt file per fasta file analyzed. Each txt file can be read as a vector. Thus the number of lines gives the dimension of the vector.
#' @seealso envnj(), fastaconc()
#' @examples \dontrun{envfastacpp("./data_t/list.txt", 10, ".", "./MyResults")}
#' @importFrom utils read.csv2
#' @export

envfascpp <- function(path, r = 10, exefile, outfile){

  l <- read.csv2(paste(path, "/list.txt", sep = ""), header = FALSE)
  l <- unlist(l)
  for (i in 1:length(l)){
    t <- strsplit(l[i], split = "/")[[1]]
    t <- t[length(t)]
    sp <- strsplit(t, "\\.")[[1]][1]
    print(paste("## --- ", i, " --- ", sp, " --- ##", sep = ""))

    com <- paste(exefile, "/envfas ", path, "/", l[i], " ", r,
                 " ", outfile, "/v_", sp, ".txt", sep = "")

    system(com)
  }
}

## ---------------------------------------------------- ##
##      vect2tree(path, metric, clustering)        ##
## ---------------------------------------------------- ##
#' Convert a Set of Vectors into a Tree
#' @description Converts a set of vectors into a tree.
#' @usage vect2tree(path, metric = "cosine", clustering = "nj")
#' @param path path to the working directory. This directory must contain a txt file per vector and an additional txt file named vlist.txt that provides the names (one per line) of the  vector txt files.
#' @param metric character string indicating the metric (see metrics() to see the methods allowed).
#' @param clustering string indicating the clustering method, either "nj" or "upgma".
#' @details This function computes the distance matrix and builds the corresponding tree.
#' @return a list with two elements: a distance matrix and a tree.
#' @seealso envnj(), fastaconc(), envfascpp()
#' @examples \dontrun{vec2tree("./data_t")}
#' @importFrom ape nj
#' @importFrom ape plot.phylo
#' @importFrom phangorn upgma
#' @export

vect2tree <- function(path, metric = "cosine", clustering = "nj"){

  oldwd <- getwd()
  on.exit(setwd(oldwd))
  setwd(path)

  l <- read.csv2("vlist.txt", header = FALSE)
  l <- unlist(l)

  # First vector/species to initialize the dataframe
  t <- strsplit(l[1], split = "/")[[1]]
  t <- t[length(t)]
  sp <- strsplit(t, "\\.")[[1]][1]
  print(paste("## --- ", 1, " --- ", sp, " --- ##", sep = ""))
  v <- read.csv2(l[1], header = FALSE)
  names(v) <- sp

  # Remaining vectors/species
  for (i in 2:length(l)){
    t <- strsplit(l[i], split = "/")[[1]]
    t <- t[length(t)]
    sp <- strsplit(t, "\\.")[[1]][1]
    print(paste("## --- ", i, " --- ", sp, " --- ##", sep = ""))
    u <- read.csv2(l[i], header = FALSE)
    names(u) <- sp
    v <- cbind(v,u)
  }
  rm(u)

  ## Distance matrix and Tree
  d <- metrics(vset = v, method = metric)
  if (clustering == 'nj'){
    t <- ape::nj(d)
  } else if (clustering == 'upgma'){
    t <- phangorn::upgma(d)
  } else {
    warning("NJ method has been selected")
  }

  ape::plot.phylo(t, use.edge.length = FALSE, cex = 0.5)

   return(list(d, t))
}




