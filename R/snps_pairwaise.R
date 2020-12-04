#' Direct SNPs pairwaise.
#'
#' This function calculates All-vs-All SNPs number. At different of \emph{core_snps_matrix} this function
#' calculates the raw number of SNPs among each pair of sequences avoiding the bias of the core genome size and/or
#' any reference.
#' 
#' @param file_list 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' @import foreachs
#' @import doParallel
#' @import dplyr
#' @import tidyr
#' @import stringr
#' 
snps_pairwaise <- function(file_list)
{
  n_cores = detectCores()-1
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  file_list <- file_list %>% as_tibble() %>% rename(File = 1)
  
  #results = data.frame()
  
  results <- foreach (i = file_list$File, .combine = "rbind", .packages = "stringr") %:%
    foreach(j = file_list$File, .combine = "rbind",.packages = "stringr") %dopar%{
      
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
  
  stopCluster(cl)
  results %>% 
    as_tibble() %>% 
    mutate(totalSNPs = as.numeric(totalSNPs)) %>% 
    mutate(totalGSNPs = as.numeric(totalGSNPs)) %>%
    mutate(AlignmentLength = as.numeric(AlignmentLength)) %>% return()
}

