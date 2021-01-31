#' Core SNPs matrix
#'
#' This function calculates the number of SNPs among samples in core genome.
#'
#' @param data A core_genome object.
#' @param norm Number of SNPs normalize by the length of the alingment?
#' @param rm.gaps remove gaps from the alignment?
#'
#' @return A square matrix
#' @export
#'
#' 
#' @import stringdist
#' @import magrittr
#'
#'


core_snps_matrix <- function(data, norm =T, rm.gaps = F){

  if(is(data,"core_genome"))
  {

    tmp <- data$core_genome$Seq %>% unlist()

    if(rm.gaps)
    {
      cmat <- str_split(tmp, pattern = "", simplify = T)
      gap.frac <- colSums(cmat == "-")

      idx.keep <- gap.frac==0

      tmp <- apply(cmat[,idx.keep], 1, paste, collapse="")
    }

    res<- stringdistmatrix(tmp, method="hamming") %>% as.matrix()
    colnames(res)<- gsub(">","",data$core_genome$Genomes)
    rownames(res) <- gsub(">","",data$core_genome$Genomes)

    if(norm)
    {
      res <- res / (nchar(data$core_genome$Seq[[1]])/1e6)
      return(res)
    }else{
      return(res)
    }
  }else if(is(data,"core_snp_genome"))
  {
    res <- stringdistmatrix(data$alignment$Seq, method="hamming") %>% as.matrix()
    colnames(res) <- data$alignment$Genomes
    rownames(res) <- data$alignment$Genomes

    if(norm)
    {
      align_length <- data$bed %>% mutate(length = END-START)
      align_length = sum(align_length$length)

      res <- res/(align_length/1e6)
      return(res)
    }else{
      return(res)
    }
  }else{
    stop("Error. data must be a core_genome or core_snp_genome object")
  }

}
