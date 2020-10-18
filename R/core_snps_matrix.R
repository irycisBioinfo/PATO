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
#'
#'


core_snps_matrix <- function(data, norm =T){

  if(!is(data,"core_genome"))
  {
    stop("data must be a core_genome object")
  }


  res <- matrix(0L, nrow =length(data$core_genome$Genomes) , ncol = length(data$core_genome$Genomes))


  for(i in 1:length(data$core_genome$Genomes)){
    print(i)

    for(j in i:length(data$core_genome$Genomes))
    {
      res[i,j] = stringdist::stringdist(data$core_genome$Seq[[i]],data$core_genome$Seq[[j]], method ="hamming")
      res[j,i] = res[i,j]
    }
  }

  colnames(res) <- gsub(">","",data$core_genome$Genomes)
  rownames(res) <- gsub(">","",data$core_genome$Genomes)

  if(norm)
  {
    res <- res / (nchar(data$core_genome$Seq[[1]])/1e6)
    return(res)
  }else{
    return(res)
  }

}
