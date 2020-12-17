#' Direct SNPs pairwaise.
#'
#' This function calculates All-vs-All SNPs number. At different of \emph{core_snps_matrix} this function
#' calculates the raw number of SNPs among each pair of sequences avoiding the bias of the core genome size and/or
#' any reference.
#'
#' @param file_list
#' @param type
#' @param n_cores
#'
#' @return
#' @export
#'
#' @examples
#'
#' @import foreach
#' @import doParallel
#' @import dplyr
#' @import tidyr
#' @import stringr
#'
snps_pairwaise <- function(file_list,type,n_cores)
{
  if(is(file_list,"gff_list"))
  {
    if(missing(type))
    {
      stop("type must be declared for gff_list objects")
    }else if(type == "nucl")
    {
      file_list = dir(paste(file_list$path,"/ffn",sep = "", collapse = ""),full.names = T) %>% as_tibble()%>% rename(File = 1)
    }else if(type =="wgs"){
      file_list = dir(paste(file_list$path,"/fna",sep = "", collapse = ""),full.names = T) %>% as_tibble()%>% rename(File = 1)
    }
  }else{
    file_list <- file_list %>% as_tibble() %>% rename(File = 1)
  }

  folderName = paste(getwd(),"/",md5(paste(file_list$File, sep = "",collapse = "")),"_spw",sep = "",collapse = "")

  if(!dir.exists(folderName))
  {
    dir.create(folderName)
  }else{
    system(paste("rm -r ",folderName))
    dir.create(folderName)
  }

  if(missing(n_cores))
  {
    n_cores = detectCores()-1
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  print("Aligning genomes")
  foreach (i = file_list$File) %dopar%{

    for(j in 1:nrow(file_list))
    {
      system(paste(minimap2," -cx asm20 -t 2 --cs ",i," ",file_list$File[j]," > ",folderName,"/",basename(i)"_",basename(file_list$File[j]),".paf",collapse = "", sep = ""), ignore.stderr = T)
      system(paste("sort -k6,6 -k8,8n ",folderName,"/",basename(i),".paf > asm.srt.paf", sep = "", collapse = ""),ignore.stderr = T)
      system(paste(k8," ",paftools," call asm.srt.paf > ",folderName,"/",basename(i),".vcf",collapse = "",sep = ""),ignore.stderr = T)
      system(paste(k8," ",paftools," splice2bed ",folderName,"/",basename(i),".paf"," > ",folderName,"/",basename(i),".bed",collapse = "",sep = ""),ignore.stderr = T)
    }
  }


}

