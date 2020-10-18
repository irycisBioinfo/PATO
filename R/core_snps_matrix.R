#' Core SNPs matrix
#'
#' This function calculates the number of SNPs among samples in core genome.
#'
#' @param data A core_genome object.
#' @param norm Number of SNPs normalize by the length of the alingment?
#'
#' @return A square matrix
#' @export
#'
#' @examples
#' @import stringdist
#' @import magrittr
#'
#'


core_snps_matrix <- function(data, norm =T){

  if(!is(data,"core_genome"))
  {
    stop("data must be a core_genome object")
  }

  res = data$core_genome$Seq %>% unlist() %>% stringdistmatrix(., method="hamming") %>% as.matrix()

  if(norm)
  {
    res <- res / (nchar(data$core_genome$Seq[[1]])/1e6)
    return(res)
  }else{
    return(res)
  }

}
