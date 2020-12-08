#' Title
#'
#' @param file_list
#' @param ref
#' @param type
#'
#' @return
#' @export
#'
#' @examples
mummer_core_snps <- function(file_list,ref,type)
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

  folderName = paste(getwd(),"/",md5(paste(file_list$File, sep = "",collapse = "")),"_mummer",sep = "",collapse = "")
  if(!dir.exists(folderName))
  {
    dir.create(folderName,)
  }



  if(missing(ref))
  {
    ref = file_list$File[runif(n = 1,min = 1,max = nrow(file_list))]
  }
  print(paste("Reference: ",ref,sep = "",collapse = ""))

  n_cores = detectCores()-1

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  foreach (i = file_list$File) %dopar%{

      system2(paste("dnadiff -p ",basename(i)," ",ref," ",i, sep = "", collapse = ""))
      system2(paste("mv ",basename(i),".* ",folderName,sep = "", collapse = ""))

  }

  stopCluster(cl)

  coords_list = dir(folderName,full.names = T) %>% grep(".1coords",.,value = T) %>% as_tibble() %>% rename(files = 1)

  c = mummer_to_data.frame(coords_list$files[1])

  for(i in 2:nrow(coords_list))
  {
    print(i)
    c = semi_join(c,mummer_to_data.frame(coords_list$files[i]))
  }


  vcf_list = list.files(folderName, pattern = ".snps", full.names = T) %>% as_tibble() %>% rename(files = 1)


  snps_table = tibble()
  for(i in 1:nrow(vcf_list))
  {
    if(file.size(vcf_list$files[i]) !=0)
    {
      tmp = fread(vcf_list$files[i]) %>% as_tibble()
      colnames(tmp) = c("pos","ref","alt","pos_query","buff","dist","r_length","q_length","rep_ref","rep_query","reference","query")
      tmp$sample = basename(vcf_list$files[i])
      snps_table = bind_rows(snps_table,tmp)
    }
  }

  snps_matrix = snps_table %>%
    filter(ref != ".") %>%
    filter(alt != ".") %>%
    select(reference,pos,ref,alt,sample) %>%
    spread(sample,alt) %>% semi_join(cov, by =c("pos" = "pos", "reference" = "contig"))

  for(i in 1:nrow(snps_matrix))
  {
    for(j in 1:ncol(snps_matrix))
    {
      snps_matrix[i,j] = ifelse(is.na(snps_matrix[i,j]), snps_matrix$ref[i], snps_matrix[i,j])
    }

  }

  return(snps_matrix)
}
