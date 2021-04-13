#' Direct SNPs pairwaise.
#'
#' This function calculates All-vs-All SNPs number. At different of \emph{core_snps_matrix} this function
#' calculates the raw number of SNPs among each pair of sequences avoiding the bias of the core genome size and/or
#' reference.
#' However, this function is an \emph{O(N²)} so can be very slow for large datasets.
#'
#' @param file_list Data frame with the full path to the nucleotide genome files (gene or genomes) or a \emph{gff_list} object.
#' @param type Just for \emph{gff_list} objects. You must especified if you want to use whole genome sequences "wgs" or genes "nucl"
#' @param n_cores Number of cores to use.
#' @param norm The output is normalized by the length of the alignment. If \emph{norm} is TRUE then the output in "SNPs per Megabase".
#' If \emph{norm} is FALSE then thr output is the raw number of SNPs.
#'
#' @return
#' @export
#'
#'
#'
#' @import foreach
#' @import doParallel
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @importFrom data.table fread
#'
snps_pairwaise <- function(file_list,type,n_cores, norm =T)
{

  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    minimap2 = system.file("minimap2",package = "pato")
    k8 = system.file("k8",package = "pato")
    paftools = system.file("paftools.js", package="pato")

  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS
    minimap2 = system.file("minimap2",package = "pato")
    k8 = system.file("k8",package = "pato")
    paftools = system.file("paftools.js", package="pato")
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
    n_cores = detectCores()/2
  }else{
    n_cores = round(n_cores/2)
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  print("Aligning genomes")
  foreach (i = 1:nrow(file_list)) %dopar%
  {
    ref = file_list$File[i]

    for(j in i:nrow(file_list))
    {
      query = file_list$File[j]

      output = paste(folderName,"/",basename(ref),"==",basename(query),".paf", sep = "", collapse = "")

      system(paste(minimap2," -cx asm20 -t 20 --cs ",ref," ",query," > ",output,collapse = "", sep = ""), ignore.stderr = T)
      system(paste("sort -k6,6 -k8,8n ",output," > ",output,".srt", sep = "", collapse = ""),ignore.stderr = T)
      system(paste(k8," ",paftools," call ",output,".srt > ",output,".vcf",collapse = "",sep = ""),ignore.stderr = T)  ### añadir parametros calling
      system(paste(k8," ",paftools," splice2bed ",output," > ",output,".bed",collapse = "",sep = ""),ignore.stderr = T)
    }
  }


  table_snps <- foreach(i = dir(folderName, pattern = ".vcf", full.names = T), .combine = "rbind") %dopar%
  {
    tmp <- data.table::fread(cmd = paste("grep '^V' ",i, sep = "", collapse = "")) %>% as_tibble()
    tmp$File <-basename(i)
    colnames(tmp) <- c("type","CHROM","POS","end","query_depth","mapping_quality","REF","ALT","query_name","query_start","query_end","query_orientation","File")
    tmp
  }

  table_snps =table_snps%>%
    filter(REF != "-") %>%
    filter(ALT != "-") %>%
    filter(str_length(REF) ==1) %>%
    filter(str_length(ALT) ==1) %>%
    group_by(File) %>%
    summarise(N_snps = n()) %>%
    separate(File,c("Source","Target"), sep = "==") %>% mutate(Target = gsub(".paf.vcf","",Target))



  if(norm)
  {
    table_beds <- foreach(i = dir(folderName, pattern = ".bed", full.names = T), .combine = "rbind") %dopar%
      {
        tmp <- data.table::fread(i) %>% as_tibble()
        tmp$File <-basename(i)
        colnames(tmp) <- c("chrom","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts","File")
        tmp
      }

    table_beds = table_beds %>%
      separate(File,c("Source","Target"), sep = "==") %>%
      mutate(Target = gsub(".paf.bed","",Target)) %>%
      mutate(SizeSource = end-start) %>%
      group_by(Source,Target) %>%
      summarise(SizeSource = sum(SizeSource)) %>% ungroup()

    table_norm = inner_join(table_snps,table_beds) %>% mutate(Norm = round(N_snps*1e6/SizeSource)) %>% select(Source,Target,Norm)
    table_norm2 = table_norm %>%
      rename(Source = Target, Target = Source)
    bind_rows(table_norm,table_norm2) %>%
      spread(Target,Norm,fill = 0) %>%
      column_to_rownames("Source")%>%
      as.matrix() %>%
      return()




  }else{
    table_snps2 = table_snps %>% rename(Source = Target, Target = Source)
    bind_rows(table_snps,table_snps2) %>%
      select(Source,Target,N_snps) %>%
      spread(Target,N_snps,fill = 0) %>%
      column_to_rownames("Source") %>%
      as.matrix() %>%
      return()

  }
}

