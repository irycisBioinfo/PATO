
#' Screening genes markers
#'
#' This function annotate virulence factors, antibiotic resitance genes and/or biocide and metals
#' of the genomes (\code{files}). The annotation is performed using
#' \strong{mmseqs2} software for protein sequences (\link{https://github.com/soedinglab/MMseqs2})
#' or \strong{Minimap2} (\link{https://github.com/lh3/minimap2}) for nucleotide (wgs or nucl).
#'
#' Databases
#' + \strong{Virulence Factor DataBase (Set_A and Set_B)} (\link{http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi})
#' + \strong{ResFinder} (\link{https://cge.cbs.dtu.dk/services/ResFinder/}).
#' + \strong{BacMet} (\link{http://bacmet.biomedicine.gu.se/}).
#'
#' The function can re-use the
#' previous computational steps of \code{mmseqs} or create a new index database from
#' the files. Re-use option shorten the computational time. This method use the algorithm
#' \strong{search} of \strong{mmseqs2} so it olny return high identity matchs.
#'
#'
#' @param data A \code{mmseq} object
#' @param type user must be specified if the data set is nucleotide or protein.
#' @param database A vector with the query databases:
#' \itemize{
#' \item \emph{AbR}: Antibiotic resistance database (ResFinder)
#' \item \emph{VF_A} VFDB core dataset (genes associated with experimentally verified VFs only)
#' \item \emph{VF_B} VFDB full dataset (all genes related to known and predicted VFs)
#' \item \emph{bacmet} BacMet dataset (genes associated with Biocide and Metal resistance)
#' }
#' @param query "all" or "accessory". It perform the annotation from whole
#' protein dataset or just from the accessory
#'
#' @return A \code{data.frame} with the annotation information\cr
#' \itemize{
#'  \item \emph{Genome}: Genome query
#'  \item \emph{Protein}: Proteins query
#'  \item \emph{target}: Protein subject (AbR o VF)
#'  \item \emph{pident}: Percentage of identical matches
#'  \item \emph{alnlen}: Alingment length
#'  \item \emph{mismatch}: number of mismatchs
#'  \item \emph{gapopen}: number of gaps
#'  \item \emph{qstart}: query start alingment
#'  \item \emph{qend}: query end alingment
#'  \item \emph{tstart}: target start alingment
#'  \item \emph{tend}:  target end alingment
#'  \item \emph{evalue}: evalue
#'  \item \emph{bits}: Bitscore
#'  \item \emph{DataBase}: Database (AbR, VF_A or VF_B, bacmet)
#'  \item \emph{Gene}: Gene name
#'  \item \emph{Description}: Functional annotation (VF) or category (AbR)
#' }
#' @note Keep in mind that the results from \emph{accesory} are based on the
#' annotation of the representative protein of the homologous cluster and therefore
#' does not mean that all the genomes have the same allele of the gene.
#'
#' @export
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import doParallel
#'
screening <- function(data, type = "nucl", database =c("AbR","VF_A","VF_B"), query = "all", n_cores)
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

  if(!(is(data,"mmseq") | is(data,"gff_list")))
  {
    stop("Error data must be a mmseq or gff_list object")
  }



  if(!dir.exists(data$path))
  {
    stop("There is an error with the input data. Please be sure that you are in the right path or re-do de mmseq object")
  }

  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    proc_cpu = readLines("/proc/cpuinfo")

    if(sum(grep("avx2",proc_cpu,ignore.case = TRUE)))
    {
      mmseqPath = system.file("mmseqs.avx2", package = "pato")
    }else{
      mmseqPath = system.file("mmseqs.sse41", package = "pato")
    }
  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS


    if(grepl("AVX2",system("sysctl -a | grep 'AVX2'", intern = T)))
    {
      mmseqPath = system.file("mmseqs.macos.avx2", package = "pato")
    }else{
      mmseqPath = system.file("mmseqs.macos.sse41", package = "pato")
    }
  }else{
    stop("Error, OS not supported.")
  }

  if(missing(n_cores))
  {
    n_cores = detectCores()
  }

  if(type == "prot")
  {
    resfinder_path <- system.file("DB/resfinder_prot", package = "pato")
    vf_A_path <- system.file("DB/VFDB_setA_prot", package = "pato")
    vf_B_path <- system.file("DB/VFDB_setB_prot", package = "pato")
    bacmet <- system.file("DB/bacmet", package = "pato")
    annot <- read.delim(system.file("DB/annot.data", package = "pato"), stringsAsFactors = FALSE, header = TRUE, sep = "\t")

  }else if(type =="nucl")
  {
    resfinder_path <- system.file("DB/resfinder_nucl", package = "pato")
    vf_A_path <- system.file("DB/VFDB_setA_nucl", package = "pato")
    vf_B_path <- system.file("DB/VFDB_setB_nucl", package = "pato")
    bacmet <- system.file("DB/bacmet", package = "pato")
    annot <- read.delim(system.file("DB/annot.data", package = "pato"), stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  }else if(type =="wgs")
  {
    resfinder_path <- system.file("DB/resfinder_nucl", package = "pato")
    vf_A_path <- system.file("DB/VFDB_setA_nucl", package = "pato")
    vf_B_path <- system.file("DB/VFDB_setB_nucl", package = "pato")
    annot <- read.delim(system.file("DB/annot.data", package = "pato"), stringsAsFactors = FALSE, header = TRUE, sep = "\t")

  }else{
    stop("Error in data type selection: please specify 'nucl', 'wgs' or 'prot'")
  }

  origin <- getwd()
  setwd(data$path)
  on.exit(setwd(origin))

  if(type == "wgs")
  {
    query ="wgs"
  }

  if(query =="all")
  {
    if(type == "prot")
    {
      rep = "all.mmseq"
    }else if (type == "nucl")
    {
      rep = "all.rnm"
    }

  }else if (query == "accessory")
  {
    if(type == "prot")
    {
      print(paste(mmseqPath," createdb all.representatives.fasta all.representatives.mm",sep = "",collapse = "")) %>% system()
      print(paste(mmseqPath," createindex all.representatives.mm tmpDir",sep = "",collapse = "")) %>% system()
      rep = "all.representatives.mm"
    }else if (type =="nucl")
    {
      rep = "all.representatives.fasta"
    }
  }else if (query == "wgs")
  {

    if(file.exists("commands.txt"))
    {
      file.remove("commands.txt")
    }

    files_wgs <- dir(paste0(data$path,"/fna/"), pattern = "fna",full.names = T) %>% as.data.frame()

    if(!file.exists(paste(data$path,"/all.rnm",sep = "",collapse = "")))
    {

      for (i in files_wgs[,1])
      {
        if(grepl("gz",i[1]))
        {
          write(paste("zcat ",i," | perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' >> ",data$path,"/all.rnm \n", collapse = "",sep = ""),
                file = "commands.txt",
                append = T)
        }else{
          write(paste("perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' ",i," >> ",data$path,"/all.rnm \n", collapse = "",sep = ""),
                file = "commands.txt",
                append = T)
        }
      }
      system(paste(Sys.getenv("SHELL")," commands.txt",collapse = "",sep = ""))
    }
    rep = paste(data$path,"/all.rnm",sep = "",collapse = "")

  }else{
    stop("Error: query must be 'all' or 'accessory'")
  }

  results <- data.frame()

  if("AbR" %in% database)
  {
    if(file.exists("abr.out.1"))
    {
      system("rm abr.out*")
    }
    if(type =="nucl" | type =="wgs")
    {
      #minimap2 -t 20 -c --cs -x asm10 ~/R/x86_64-pc-linux-gnu-library/4.0/pato/DB/resfinder_nucl all.rnm > all.paf
      print(paste(minimap2," -c --cs -x asm10 -t ",n_cores," ",resfinder_path," ",rep," > abr.paf", sep = "", collapse = "")) %>% system()

      Sys.sleep(1)

      tmp <- read_table("abr.paf", col_names = F)
      tmp <- tmp %>% select(X1:X19)
      colnames(tmp) <- c("qname","qlength","qstart","qend","strand",
                         "target","tlength","tstart","tend","nmatches",
                         "alnlength","mapq","NM","ms","AS","nn","tp","cm","s1")
      tmp <- tmp %>%
        mutate(NM = str_replace(NM,"NM:i:","") %>% as.numeric()) %>%
        mutate(ms = str_replace(ms,"ms:i:","") %>% as.numeric()) %>%
        mutate(AS = str_replace(AS,"AS:i:","") %>% as.numeric()) %>%
        mutate(nn = str_replace(nn,"nn:i:","") %>% as.numeric()) %>%
        mutate(tp = str_replace(tp,"tp:A:","")) %>%
        mutate(cm = str_replace(cm,"cm:i:","") %>% as.numeric()) %>%
        mutate(s1 = str_replace(s1,"s1:i:","") %>% as.numeric()) %>%
        mutate(identity = nmatches/alnlength) %>%
        mutate(coverage = (tend - tstart)/tlength) %>%
        rename(score = AS) %>%
        rename(mismatches = NM) %>%
        mutate(aln = cumsum(ifelse(tp=="P",1,0))) %>%
        group_by(aln) %>%
        slice_max(order_by = score,with_ties = T) %>%
        ungroup() %>%
        group_by(qname,target,aln) %>%
        summarise_all(first)

      results <- bind_rows(results,tmp)


    }else if (type =="prot")
    {
      print(paste(mmseqPath," map ",rep," ",resfinder_path," abr.out tmpDir", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)

      print(paste(mmseqPath," convertalis ",rep," ",resfinder_path," abr.out abr.tsv", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
      tmp<- read.table("abr.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
      Sys.sleep(1)
      colnames(tmp) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")
      tmp <- tmp %>% mutate(target = gsub("_\\d+$","",target))
      results <- bind_rows(results,tmp)
    }
  }
  if("VF_A" %in% database)
  {
    if(file.exists("vf_a.out.1"))
    {
      system("rm vf_a.out*")
    }

    if(type =="nucl" | type =="wgs")
    {
      print(paste(minimap2," -c --cs -x asm10 -t ",n_cores," ",vf_A_path," ",rep," > vfa.paf", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)

      tmp <- read_table("vfa.paf", col_names = F)
      tmp <- tmp %>% select(X1:X19)
      colnames(tmp) <- c("qname","qlength","qstart","qend","strand",
                         "target","tlength","tstart","tend","nmatches",
                         "alnlength","mapq","NM","ms","AS","nn","tp","cm","s1")
      tmp <- tmp %>%
        mutate(NM = str_replace(NM,"NM:i:","") %>% as.numeric()) %>%
        mutate(ms = str_replace(ms,"ms:i:","") %>% as.numeric()) %>%
        mutate(AS = str_replace(AS,"AS:i:","") %>% as.numeric()) %>%
        mutate(nn = str_replace(nn,"nn:i:","") %>% as.numeric()) %>%
        mutate(tp = str_replace(tp,"tp:A:","")) %>%
        mutate(cm = str_replace(cm,"cm:i:","") %>% as.numeric()) %>%
        mutate(s1 = str_replace(s1,"s1:i:","") %>% as.numeric()) %>%
        mutate(identity = nmatches/alnlength) %>%
        mutate(coverage = (tend - tstart)/tlength) %>%
        rename(score = AS) %>%
        rename(mismatches = NM) %>%
        mutate(aln = cumsum(ifelse(tp=="P",1,0))) %>%
        group_by(aln) %>%
        slice_max(order_by = score,with_ties = T) %>%
        ungroup() %>%
        group_by(qname,target,aln) %>%
        summarise_all(first)

      results <- bind_rows(results,tmp)

    }else{
      print(paste(mmseqPath," map ",rep," ",vf_A_path," vf_a.out tmpDir", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
      print(paste(mmseqPath," convertalis ",rep," ",vf_A_path," vf_a.out vf_a.tsv", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
      tmp <- read.table("vf_a.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
      Sys.sleep(1)
      colnames(tmp) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")

      Sys.sleep(1)
      results <- bind_rows(results,tmp)
    }



  }
  if("VF_B" %in% database)
  {
    if(file.exists("vf_b.out.1"))
    {
      system("rm vf_b.out*")
    }

    if(type =="nucl" | type =="wgs")
    {
      print(paste(minimap2," -c --cs -x asm10 -t ",n_cores," ",vf_B_path," ",rep," > vfa.paf", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)

      tmp <- read_table("vfa.paf", col_names = F)
      tmp <- tmp %>% select(X1:X19)
      colnames(tmp) <- c("qname","qlength","qstart","qend","strand",
                         "target","tlength","tstart","tend","nmatches",
                         "alnlength","mapq","NM","ms","AS","nn","tp","cm","s1")
      tmp <- tmp %>%
        mutate(NM = str_replace(NM,"NM:i:","") %>% as.numeric()) %>%
        mutate(ms = str_replace(ms,"ms:i:","") %>% as.numeric()) %>%
        mutate(AS = str_replace(AS,"AS:i:","") %>% as.numeric()) %>%
        mutate(nn = str_replace(nn,"nn:i:","") %>% as.numeric()) %>%
        mutate(tp = str_replace(tp,"tp:A:","")) %>%
        mutate(cm = str_replace(cm,"cm:i:","") %>% as.numeric()) %>%
        mutate(s1 = str_replace(s1,"s1:i:","") %>% as.numeric()) %>%
        mutate(identity = nmatches/alnlength) %>%
        mutate(coverage = (tend - tstart)/tlength) %>%
        rename(score = AS) %>%
        rename(mismatches = NM) %>%
        mutate(aln = cumsum(ifelse(tp=="P",1,0))) %>%
        group_by(aln) %>%
        slice_max(order_by = score,with_ties = T) %>%
        ungroup() %>%
        group_by(qname,target,aln) %>%
        summarise_all(first)

      results <- bind_rows(results,tmp)

    }else{
      print(paste(mmseqPath," map ",rep," ",vf_B_path," vf_b.out tmpDir", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
      print(paste(mmseqPath," convertalis ",rep," ",vf_B_path," vf_b.out vf_b.tsv", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
      tmp <- read.table("vf_b.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
      Sys.sleep(1)
      colnames(tmp) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")

      Sys.sleep(1)
      results <- bind_rows(results,tmp)
    }

  }
  if("bacmet" %in% database)
  {
    if(file.exists("bacmet.out.1"))
    {
      system("rm bacmet.out*")
    }
    if(type =="nucl")
    {
      print(paste(mmseqPath," search ",rep," ",bacmet," bacmet.out tmpDir --search-type 3 ", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
    }else{
      print(paste(mmseqPath," map ",rep," ",bacmet," bacmet.out tmpDir", sep = "", collapse = "")) %>% system()
      Sys.sleep(1)
    }


    print(paste(mmseqPath," convertalis ",rep," ",bacmet," bacmet.out bacmet.tsv", sep = "", collapse = "")) %>% system()
    Sys.sleep(1)
    tmp <- read.table("bacmet.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
    Sys.sleep(1)
    colnames(tmp)<- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")
    Sys.sleep(1)
    results <- bind_rows(results,tmp)
  }


  if(type == "prot")
  {
    results <- inner_join(results,annot, by = c("target" = "ID")) %>% separate(query,c("Genome","Protein"), sep = "#")
    if(!("VF_B" %in% database)){
      results <- results %>% filter(DataBase != "VFB")
    }
  }else{
    results <- inner_join(results,annot, by = c("target" = "ID")) %>% separate(qname,c("Genome","Protein"), sep = "#")
    if(!("VF_B" %in% database)){
      results <- results %>% filter(DataBase != "VFB")
    }
  }

  return(results)
}



