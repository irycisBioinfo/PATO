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

  if(file.exists(file))
  {
    file.remove(file)
  }

  if(is(core_data,"core_genome"))
  {
    for(i in 1:length(core_data$core_genome$Genomes))
    {
      write(core_data$core_genome$Genomes[i],file,append = TRUE)
      write(core_data$core_genome$Seq[i],file, append = TRUE)
    }
  }else if(is(core_data,"core_snp_genome")){
    for(i in 1:nrow(core_data$alignment))
    {
      write(paste(">",core_data$alignment$Genomes[i],sep = "",collapse = ""),file,append = TRUE)
      write(core_data$alignment$Seq[i],file, append = TRUE)

    }
  }else{
    stop("core_data must be a 'core_genome' object")
  }
}

