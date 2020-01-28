#' Export core-genome alignment to Multi-alignment FASTA file.
#'
#' @param core_data A core_genome object
#' @param file Path for alignment ouput file
#'
#' @export
#'
#' @examples
export_core_to_fasta <- function(core_data,file)
{
  if(is(core_data,"core_genome"))
  {
    for(i in 1:length(core_data$Genomes))
    {
      write(core_data$Genomes[i],file,append = TRUE)
      write(core_data$Seq[i],file, append = TRUE)
    }
  }else{
    stop("core_data")
  }
}

