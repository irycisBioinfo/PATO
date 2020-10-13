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
#' #' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#' @import stringr
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#'


core_snps_matrix <- function(data, norm =T){

  if(!is(data,"core_genome"))
  {
    stop("data must be a core_genome object")
  }

  # n_cores = detectCores()
  # cl <- makeCluster(n_cores)
  # registerDoParallel(cl)

  res <- matrix(0L, nrow =length(data$core_genome$Genomes) , ncol = length(data$core_genome$Genomes))


  for(i in 1:length(data$core_genome$Genomes)){
    print(i)
    tmp1 <- data$core_genome$Seq[i] %>% str_split("") %>% as.matrix()
    for(j in i:length(data$core_genome$Genomes))
    {

      tmp2 <- data$core_genome$Seq[j] %>% str_split("") %>% as.matrix()
      res[i,j] = sum(tmp1[[1]] != tmp2[[1]])
      #res[i,j] = adist(data$core_genome$Seq[i],data$core_genome$Seq[j])
      res[j,i] = res[i,j]
    }
  }
  # stopCluster(cl)

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
