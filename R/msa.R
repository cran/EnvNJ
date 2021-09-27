## ------------- msa.R ------------- ##
#                                     #
#    msa.merge                        #
#    msa.tree                         #
#                                     #
## --------------------------------- ##


## ----------------------------------------------------------- ##
#                 msa.merge(data, outfile)                      #
## ----------------------------------------------------------- ##
#' Carry Out a MSA of a Set of Different Orthologous Proteins
#' @description Carries out a MSA of a set of different orthologous proteins in different species.
#' @usage  msa.merge(data, outfile = 'any')
#' @param data input data must be a dataframe where each row corresponds to a protein and each column to a species.
#' @param outfile path to the place where a fasta file is going to be saved. If 'any', no file is saved.
#' @details The input data has the same format that the input data used for EnvNJ or SVD-n-Gram methods. Thus, the name of columns must correspond to that of species.
#' @return A dataframe containing the MSA (species x position).
#' @seealso msa.tree
#' @examples \dontrun{
#' data(bovids)
#' msa.merge(bovids)}
#' @importFrom seqinr write.fasta
#' @importFrom bio3d seqbind
#' @importFrom bio3d seqaln
#' @export

msa.merge <- function(data, outfile = 'any'){

  sp <- names(data)

  for (i in 1:nrow(data)){ # for each protein
    s <- unlist(data[i, ])
    if (length(s) < 3){
      stop("Too few orthologous sequences!")
    }
    sb <- bio3d::seqbind(s[1], s[2])
    for (j in 3:length(s)){
      sb <- seqbind(sb, s[j])
    }
    a <- seqaln(sb)

    t <- a$ali
    if (i == 1){
      msa <- t
    } else {
      msa <- cbind(msa, t)
    }
  }
  if (outfile != 'any'){
    sequences <- apply(msa, 1, function(x) paste(x, collapse = ""))
    seqinr::write.fasta(sequences = as.list(sequences), names = sp, file.out = outfile)
  }

  if (file.exists("aln.fa")){
    file.remove("aln.fa")
  }
  rownames(msa) <- sp
  return(as.data.frame(msa))
}


## ----------------------------------------------------------- ##
#                 msa.tree(data, outgroup)                      #
## ----------------------------------------------------------- ##
#' Infer a tree based on a MSA
#' @description Infers a tree base on a MSA.
#' @usage  msa.tree(data, outgroup = 'any')
#' @param data input data must be a dataframe where each row corresponds to a protein and each column to a species.
#' @param outgroup when a rooted tree is desired, indicate the species to be used as outgroup.
#' @details The input data has the same format that the input data used for EnvNJ or SVD-n-Gram methods. Thus, the name of columns must correspond to that of species.
#' @return A list containing the (i) MSA, (ii) the distance matrix and (iii) the tree.
#' @seealso msa.merge
#' @examples \dontrun{
#' data(bovids)
#' msa.tree(bovids)}
#' @importFrom ape dist.aa
#' @importFrom ape nj
#' @importFrom ape root
#' @importFrom ape plot.phylo
#' @export

msa.tree <- function(data, outgroup = 'any'){
  msa <- msa.merge(data)
  d_msa <- as.matrix(ape::dist.aa(msa))
  t_msa <- ape::nj(d_msa)
  if (outgroup[1] != 'any'){
    t_msa <- ape::root(t_msa, outgroup = outgroup)
  }

  ape::plot.phylo(t_msa, cex = 0.5, use.edge.length = FALSE)

  return(list(msa, d_msa, t_msa))
}


