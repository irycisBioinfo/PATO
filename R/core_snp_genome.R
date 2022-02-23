#' Core SNP Genome
#'
#' This function find the core snp genome. Citing \emph{Torsten Seemann} :
#' \emph{If you call SNPs for multiple isolates from the same reference,
#' you can produce an alignment of "core SNPs" which can be used to build a
#' high-resolution phylogeny (ignoring possible recombination). A "core
#' site" is a genomic position that is present in all the samples. A core
#' site can have the same nucleotide in every sample ("monomorphic") or some
#' samples can be different ("polymorphic" or "variant"). If we ignore the
#' complications of "ins", "del" variant types, and just use variant sites,
#' these are the "core SNP genome".}
#'
#' This function uses \emph{minimap2} to align all the genomes to a reference genome.
#' If reference genome is not specified then \emph{core_snp_genome()} takes one
#' ramdomly. Once we have all the genomes aligned we look for the conserved regions
#' between all the genomes. Then, we call for the variants using the tool provided
#' with \emph{minimap2}: \emph{paftools}.
#'
#' Finally the SNPs al filtered by the common regions to produce the final core SNP genome.
#'
#' @param file_list Data frame with the full path to the nucleotide genome files (gene or genomes) or a \emph{gff_list} object.
#' @param n_cores Number of cores to use.
#' @param ref Reference genome (if missing, one is selected randomly)
#' @param type Just for \emph{gff_list} objects. You must especified if you want to use whole genome sequences "wgs" or genes "nucl"
#' @param x minimap preset (see details)
#' @param min_call_length min alignment length to call variants and compute coverage (expert parameters)
#' @param min_call_qual min mapping quality (expert parameters)
#' @param asm Minimap2 preseting options (see details)
#'
#'
#' @details
#' Minimap has some preset setting to map different kind of sequences.
#' - map-pb/map-ont: PacBio/Nanopore vs reference mapping
#' - ava-pb/ava-ont: PacBio/Nanopore read overlap
#' - asm5/asm10/asm20: asm-to-ref mapping, for ~0.1/1/5% sequence divergence
#' - splice: long-read spliced alignment
#' - sr: genomic short-read mapping
#'
#' We recommend to use *map-bp* (maximum number of SNPs ~ less accuracy) or
#' *asm5* (lowest SNP max accuracy), *asm10* (medium SNP, medium accuracy) or
#' *asm20* (high SNP, low accuracy)
#'
#' The rest of the presets are designed for other purposes
#'
#' @references Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191
#'
#' @return core_snp_genome object
#' @export
#'
#'
core_snp_genome <- function(file_list, n_cores, ref, type, x,min_call_length, min_call_qual, asm)
{
  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    minimap2 <- system.file("minimap2",package = "pato")
    k8 <- system.file("k8",package = "pato")
    paftools <- system.file("paftools.js", package="pato")
    bedtools <- system.file("bedtools", package="pato")

  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS
    minimap2 <- system.file("minimap2",package = "pato")
    k8 <- system.file("k8",package = "pato")
    paftools <- system.file("paftools.js", package="pato")
    bedtools <- system.file("bedtools", package="pato")

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

  if(missing(n_cores))
  {
    n_cores <- detectCores()/2
  }else{
    n_cores <- n_cores/2
  }

  if(missing(min_call_length))
  {
    min_call_length = 1000
  }

  if(missing(min_call_qual))
  {
    min_call_qual = 5
  }
  if(missing(asm))
  {
    asm = "map-pb"
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  options(warn = -1)


  print("Indexing reference genome")
  system(paste(minimap2," -d ref.mmi ",ref, sep = "", collapse = ""))
  file.copy(from = ref,to = paste0(folderName,"/reference.fasta"))

  print("Aligning genomes")
  foreach (i = file_list$File) %dopar%{

    system(paste(minimap2," -cx ",asm," -t 2 --cs=long ref.mmi ",i," > ",folderName,"/",basename(i),".paf",collapse = "", sep = ""), ignore.stderr = T)
    #system(paste(minimap2," -c -t 2 --cs=long ref.mmi ",i," > ",folderName,"/",basename(i),".paf",collapse = "", sep = ""), ignore.stderr = T)
    system(paste("sort -k6,6 -k8,8n ",folderName,"/",basename(i),".paf > ",folderName,"/",basename(i),".tmp.paf", sep = "", collapse = ""),ignore.stderr = T)
    system(paste0(k8," ",paftools," call -L ",min_call_length," -l ",min_call_length," -q ",min_call_qual," ",folderName,"/",basename(i),".tmp.paf > ",folderName,"/",basename(i),".vcf",collapse = "",sep = ""),ignore.stderr = T)
    system(paste(k8," ",paftools," splice2bed ",folderName,"/",basename(i),".paf"," > ",folderName,"/",basename(i),".tmp",collapse = "",sep = ""),ignore.stderr = T)
    system(paste0(bedtools," sort -i ",folderName,"/",basename(i),".tmp > ",folderName,"/",basename(i),".bed"))
  }
  stopCluster(cl)

  options(warn = 0)

  vcfs <- dir(folderName, pattern = ".vcf", full.names = T)
  beds <- dir(folderName, pattern = ".bed", full.names = T)


  # system(
  #    paste(
  #      "bedtools intersect -a ",
  #      beds[1],
  #      " -b ",
  #      paste(beds[-1], sep = " ", collapse = " "),
  #      " | sort -k1,1 -k2,2n | bedtools merge > final.bed",
  #      sep = "",
  #      collapse = ""
  #    )
  #  )
  system(paste0("cp ",beds[1]," ",folderName,"/final.bed"))
  system(paste0("cp ",beds[1]," ",folderName,"/ini.bed"))
  for(i in beds[-1])
  {
    system(paste0("bedtools intersect -a ",folderName,"/final.bed -b ",i," | sort -k1,1 -k2,2n | bedtools merge > ",folderName,"/tmp.bed"))
    system(paste0("cp ",folderName,"/tmp.bed ",folderName,"/final.bed"))
  }
  print("Reading Bed Files")
   bed <- read.table(paste0(folderName,"/final.bed"), stringsAsFactors = F)
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
    tmp <- fread(sep = "\t",header = F, stringsAsFactors = F,cmd = paste("grep '^V' ",vcfs[i], sep = "",collapse ="" ))
    tmp$sample <- basename(vcfs[i])
    vcf_table <- bind_rows(vcf_table,tmp)
  }
  colnames(vcf_table) <- c("type","CHROM","POS","end","query_depth","mapping_quality","REF","ALT","query_name","query_start","query_end","query_orientation","Sample")
  vcf_table <- vcf_table %>% as_tibble()

  vcf_table$POS <- as.numeric(vcf_table$POS)
  positions$POS <- as.numeric(positions$POS)

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
  result = result %>%
    ungroup() %>%
    select(-CHROM, -POS) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("Genomes") %>%
    unite(Seq,-Genomes,sep = "")
    result$Genomes = gsub(".vcf","",result$Genomes)

  output <- list(alignment = result,bed = bed,path = folderName, reference = ref)
  class(output) <-  append(class(output),"core_snp_genome" )
  return(output)





}
