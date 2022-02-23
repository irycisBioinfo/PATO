#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
eggnogmapper <- function(data, n_cores)
{

  system(paste0("emapper.py --override -o ",data$path,"/eggnog -i ",data$path,"/all.representatives.fasta"," --cpu ",n_cores))
  result <- read_tsv(paste0(data$path,"/eggnog.emapper.annotations"),comment = "##")
  return(result)

}
