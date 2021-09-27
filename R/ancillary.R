## ----.--- ancillary.R ----------- ##
#                                    #
#   aa.at                            #
#   aa.comp                          #
#   aaf                              #
#                                    #
## -------------------------------- ##

## --------------------------------------------------------------- ##
#                        aa.at(at, target)                          #
## --------------------------------------------------------------- ##
#' Residue Found at the Requested Position
#' @description Returns the residue found at the requested position.
#' @usage aa.at(at, target)
#' @param at the position in the primary structure of the protein.
#' @param target a character string specifying the sequence of the protein of interest.
#' @return Returns a single character representing the residue found at the indicated position in the indicated protein.
#' @examples aa.at(2, 'MFSQQQRCVE')
#' @seealso aa.comp()
#' @importFrom bio3d get.seq
#' @export

aa.at <- function(at, target){

 if (length(target) == 1){
      target <- strsplit(target, split="")[[1]]
 }

 if (at %in% 1:length(target)){
    return(target[at])
 } else {
    message(paste(at , " isn't a valid position for this protein", sep=""))
    return(NULL)
 }
}


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
#                          aaf(data)                            #
## ----------------------------------------------------------- ##
#' Compute the Frequency of Each Amino Acid in Each Species
#' @description Computes the frequency of each amino acid in each species.
#' @usage aaf(data)
#' @param data input data must be a dataframe (see details).
#' @details Input data must be a dataframe where each row corresponds to an individual protein, and each column identifies a species. Therefore, the columns' names of this dataframe must be coherent with the names of the OTUs being analyzed.
#' @return A dataframe providing amino acid frequencies en the set of species. Rows correspond amino acids and columns to species.
#' @seealso env.sp(), otu.vector(), otu.space()
#' @examples aaf(bovids)
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

