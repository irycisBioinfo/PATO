#' EggNogg Mapper
#'
#' this a wraper function to external software EggnoggMapper
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
eggnogmapper <- function(data, n_cores, executable_path = "emapper.py")
{


  system(paste0(executable_path, "--override -o ",data$path,"/eggnog -i ",data$path,"/all.representatives.fasta"," --cpu ",n_cores))
  result <- read_tsv(paste0(data$path,"/eggnog.emapper.annotations"),comment = "##") %>% rename(query = 1)
  return(result)

}
