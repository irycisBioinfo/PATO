#' Title
#'
#' @param file_list
#' @param n_core
#' @param ref
#'
#' @return
#' @export
#'
#' @examples
core_by_mapping <- function(file_list, n_core, ref, type)
{
  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    gsaPath = system.file("GSAlign",package = "pato")
    bwtPath = system.file("bwt_index",package = "pato")

  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS
    ### Need to compile
  }else{
    stop("Error, OS not supported.")
  }

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

  folderName = paste(getwd(),"/",md5(paste(file_list$File, sep = "",collapse = "")),"_gsa",sep = "",collapse = "")

  if(!dir.exists(folderName))
  {
    dir.create(folderName,)
  }

  if(missing(ref))
  {
    ref = file_list$File[runif(n = 1,min = 1,max = nrow(file_list))]
  }
  print(paste("Reference: ",ref,sep = "",collapse = ""))

  #origin_path = getwd()
  #setwd(folderName)
  #on.exit(setwd(origin_path))

  system(paste(bwtPath," ",ref," ",folderName,"/reference.bwt", collapse = "", sep = ""))



  n_cores <- detectCores()/2

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  options(warn = -1)


  foreach (i = file_list$File) %dopar%{

    system(paste(gsaPath," -one -i ",folderName,"/reference.bwt -q ",i," -o ",folderName,"/",basename(i),".output"," -t 2", sep = "", collapse = ""),ignore.stdout = T, ignore.stderr = T)

  }
  stopCluster(cl)

  options(warn = 0)

  mafs <- dir(folderName, pattern = ".maf", full.names = T)



  for(i in 1:length(mafs))
  {
    maf_to_data.frame(mafs[i])
  }
  beds = gsub(".maf",".bed",mafs)
  system(paste("bedtools intersect -a ",beds[1]," -b ",paste(beds[-1],sep = " ",collapse = " ")," >final.bed",sep = "",collapse = ""))
  result = read.table("final.bed")
  return(result)
}
