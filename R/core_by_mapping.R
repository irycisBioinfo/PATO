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
    minimap2 = system.file("minimap2",package = "pato")
    k8 = system.file("k8",package = "pato")
    paftools = system.file("paftools.js", package="pato")

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

  folderName = paste(getwd(),"/",md5(paste(file_list$File, sep = "",collapse = "")),"_mmi",sep = "",collapse = "")

  if(!dir.exists(folderName))
  {
    dir.create(folderName)
  }else{
    system(paste("rm -r ",folderName))
    dir.create(folderName)
  }

  if(missing(ref))
  {
    ref = file_list$File[runif(n = 1,min = 1,max = nrow(file_list))]
  }
  print(paste("Reference: ",ref,sep = "",collapse = ""))




  n_cores <- detectCores()/2

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  options(warn = -1)


  print("Indexing reference genome")
  system(paste(minimap2," -d ref.mmi ",ref, sep = "", collapse = ""))

  print("Aligning genomes")
  foreach (i = file_list$File) %dopar%{


    system(paste(minimap2," -cx asm20 -t 2 --cs ref.mmi ",i," > ",folderName,"/",basename(i),".paf",collapse = "", sep = ""), ignore.stderr = T)
    system(paste(k8," ",paftools," call ",folderName,"/",basename(i),".paf"," > ",folderName,"/",basename(i),".vcf",collapse = "",sep = ""),ignore.stderr = T)
    system(paste(k8," ",paftools," splice2bed ",folderName,"/",basename(i),".paf"," > ",folderName,"/",basename(i),".bed",collapse = "",sep = ""),ignore.stderr = T)
  }
  stopCluster(cl)

  options(warn = 0)

  vcfs <- dir(folderName, pattern = ".vcf", full.names = T)
  beds <- dir(folderName, pattern = ".bed", full.names = T)

  system(
     paste(
       "bedtools intersect -a ",
       beds[1],
       " -b ",
       paste(beds[-1], sep = " ", collapse = " "),
       " | sort -k1,1 -k2,2n | bedtools merge > final.bed",
       sep = "",
       collapse = ""
     )
   )
  print("Reading Bed Files")
   bed <- read.table("final.bed")
   colnames(bed) <- c("CHROM","START","END")


   positions <- data.frame()
   for(i in 1:nrow(bed))
   {
     positions <- bind_rows(positions, data.frame(POS = seq(bed$START[i]:bed$END[i]), CHROM = bed$CHROM[i], stringsAsFactors = F))
   }

  print("Reading VCF files")
  vcfs <- vcfs[!grepl(basename(ref),vcfs)]
  vcf_table <- data.frame()
  for(i in 1:length(vcfs))
  {
    tmp <- fread(sep = "\t",header = F, stringsAsFactors = F,cmd = paste("grep 'V' ",vcfs[i], sep = "",collapse = ))
    tmp$sample = basename(vcfs[i])
    vcf_table = bind_rows(vcf_table,tmp)
  }
  colnames(vcf_table) <- c("type","CHROM","POS","end","query_depth","mapping_quality","REF","ALT","query_name","query_start","query_end","query_orientation","Sample")
  vcf_table <- vcf_table %>% as_tibble()


  result <- semi_join(vcf_table,positions)
  result <- result %>%
    filter(str_length(REF) ==1)%>%
    filter(str_length(ALT) ==1)%>%
    filter(REF != "-")%>%
    filter(ALT != "-")%>%
    select(CHROM,POS,REF,Sample,ALT) %>%
    distinct() %>%
    group_by(CHROM,REF,POS,Sample) %>%
    mutate(Hplot = n()) %>%
    filter(Hplot ==1)%>%
    select(-Hplot) %>%
    spread(Sample,ALT) %>%
    ungroup()

  for(i in 4:ncol(result))
  {
    result[[i]] <- coalesce(result[[i]],result$REF)
  }
  RefName = basename(ref)
  colnames(result) = gsub("REF",RefName,colnames(result))
  result = result %>% ungroup() %>% select(-CHROM, -POS) %>% t() %>% as.data.frame() %>% rownames_to_column("Genomes") %>% unite(Seq,-Genomes,sep = "")

  return(list(aligment = result,vcf= vcf_table,bed = bed,path = folderName))
  #return(list(vcf= vcf_table,pos = positions))



}
