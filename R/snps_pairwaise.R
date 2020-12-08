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

  if(missing(n_cores))
  {
    n_cores = detectCores()-1
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)


  results <- foreach (i = file_list$File, .combine = "rbind", .packages = "stringr") %:%
    foreach(j = file_list$file, .combine = "rbind",.packages = "stringr") %do%{



      print(paste("dnadiff -p ",basename(i),"_",basename(j)," ",i," ",j, sep = "", collapse = ""))
      system(paste("dnadiff -p ",basename(i),"_",basename(j)," ",i," ",j, sep = "", collapse = ""))

      report = readLines(paste(basename(i),"_",basename(j),".report",sep = "", collapse = ""))

      totalSNPs = str_split(grep("TotalSNPs",report, value = TRUE),pattern = "\\s+")[[1]][2]
      totalGSNPs = str_split(grep("TotalGSNPs",report, value = TRUE),pattern = "\\s+")[[1]][2]
      AlignmentLength = str_split(report[17],pattern = "\\s+")[[1]][2]

      system(paste("rm ",basename(i),"_",basename(j),".*",sep = "", collapse = ""))
      results = c(Source =basename(i),
                  Target = basename(j),
                  totalSNPs = as.numeric(totalSNPs),
                  totalGSNPs= as.numeric(totalGSNPs),
                  AlignmentLength = as.numeric(AlignmentLength))
  }

  on.exit(stopCluster(cl))

  results %>%
    as_tibble() %>%
    mutate(totalSNPs = as.numeric(totalSNPs)) %>%
    mutate(totalGSNPs = as.numeric(totalGSNPs)) %>%
    mutate(AlignmentLength = as.numeric(AlignmentLength)) %>% return()
}

