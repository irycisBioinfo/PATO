#' Annotate Virulence Factors and Antibiotic Resistance Genes
#'
#' This function annotate virulence factors and the antibiotic resitance genes
#' of the genomes (\code{files}). The annotation is performed using
#' \strong{mmseqs2} software (\link{https://github.com/soedinglab/MMseqs2}) and the
#' databases \strong{Virulence Factor DataBase}
#' (\link{http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi}) and \strong{ResFinder}
#' (\link{https://cge.cbs.dtu.dk/services/ResFinder/}). The function can re-use the
#' previous computational steps of \code{mmseqs} or create a new index database from
#' the files. Re-use option shorten the computational time. This method use the algorithm
#' \strong{map} of \strong{mmseqs2} so it olny return high identity matchs.
#'
#'
#' @param files data.frame with the absolute path to the genome files (protein fasta file).
#' @param re_use Re-use the precomputed mmseqs2 database. The index must be in
#' the working directory
#' @param type user must be specified if the data set is nucleotide or protein.
#' @param database A vector with the query databases:
#' \itemize{
#' \item \emph{AbR}: Antibiotic resistance database (ResFinder)
#' \item \emph{VF_A} VFDB core dataset (genes associated with experimentally verified VFs only)
#' \item \emph{VF_B} VFDB full dataset (all genes related to known and predicted VFs)
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
#'  \item \emph{DataBase}: Database (AbR, VF_A or VF_B)
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
#' @import dtplyr
#' @import data.table
#'
annotate <- function(files, re_use = TRUE, type = "nucl", database =c("AbR","VF_A","VF_B"), query = "all")
{

  if(sum(grep("avx2",system("cat /proc/cpuinfo",intern = TRUE),ignore.case = TRUE)))
  {
    mmseqPath = system.file("mmseqs.avx2", package = "PATO")
  }else{
    mmseqPath = system.file("mmseqs.sse41", package = "PATO")
  }
  if(type == "prot")
  {
    resfinder_path <- system.file("annotation/resfinder_prot", package = "PATO")
    vf_A_path <- system.file("annotation/VFDB_setA_prot", package = "PATO")
    vf_B_path <- system.file("annotation/VFDB_setB_prot", package = "PATO")
    annot <- read.delim(system.file("annotation/annot.data", package = "PATO"), stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  }else if(type =="nucl")
  {
    resfinder_path <- system.file("annotation/resfinder_nucl", package = "PATO")
    vf_A_path <- system.file("annotation/VFDB_setA_nucl", package = "PATO")
    vf_B_path <- system.file("annotation/VFDB_setB_nucl", package = "PATO")
    annot <- read.delim(system.file("annotation/annot.data", package = "PATO"), stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  }else{
    stop("Error in data type selection: please specify 'nucl' or 'prot'")
  }

  results <- data.frame()

  if(query =="all")
  {
    if(!re_use)
    {
      system("rm *.rnm")
      system("rm -r tmpDir")
      for (i in files[,1])
      {

        if(grepl("gz",i[1]))
        {
          print(paste("zcat ",i," | perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' >> all.rnm", collapse = "",sep = ""), quote = FALSE) %>% system()
        }else{
          print(paste("perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' ",i," >> all.rnm", collapse = "",sep = ""), quote = FALSE) %>% system()

        }

      }

      paste(mmseqPath," createdb all.rnm all.mmseq",sep = "",collapse = "") %>% system()
      paste(mmseqPath," createindex all.mmseq tmpDqir",sep = "",collapse = "") %>% system()

    }

    if(!file.exists("all.mmseq"))
    {
      stop("all.mmseq file not found. Try option: re_use = FALSE")
    }

    if(!prod(database %in% c("AbR","VF_A","VF_B")))
    {
      stop("Error: Only AbR, VF_A and VF_B are available")
    }
    rep = "all.mmseq"
  }else if (query == "accessory")
  {
    if(file.exists("all.representatives.fasta"))
    {
      paste(mmseqPath," createdb all.representatives.fasta all.representatives.mm",sep = "",collapse = "") %>% system()
      paste(mmseqPath," createindex all.representatives.mm tmpDqir",sep = "",collapse = "") %>% system()

      rep = "all.representatives.mm"
    }else{
      stop("Does not exists 'all.representatives.fasta' file in your workdirectory")
    }
  }else{
    stop("Error: query must be 'all' or 'accessory'")
  }


  if("AbR" %in% database)
  {
    if(file.exists("abr.out"))
    {
      system("rm abr.out")
    }

    paste(mmseqPath," map ",rep," ",resfinder_path," abr.out tmpDir", sep = "", collapse = "") %>% system()
    paste(mmseqPath," convertalis ",rep," ",resfinder_path," abr.out abr.tsv", sep = "", collapse = "") %>% system()
    tmp<- read.table("abr.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
    colnames(tmp) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")

    results <- bind_rows(results,tmp)
  }
  if("VF_A" %in% database)
  {
    if(file.exists("vf_a.out"))
    {
      system("rm vf_a.out")
    }
    paste(mmseqPath," map ",rep," ",vf_A_path," vf_a.out tmpDir", sep = "", collapse = "") %>% system()
    paste(mmseqPath," convertalis ",rep," ",vf_A_path," vf_a.out vf_a.tsv", sep = "", collapse = "") %>% system()
    tmp <- read.table("vf_a.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
    colnames(tmp) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")
    results <- bind_rows(results,tmp)

  }
  if("VF_B" %in% database)
  {
    if(file.exists("vf_b.out"))
    {
      system("rm vf_b.out")
    }
    paste(mmseqPath," map ",rep," ",vf_B_path," vf_b.out tmpDir", sep = "", collapse = "") %>% system()
    paste(mmseqPath," convertalis ",rep," ",vf_B_path," vf_b.out vf_b.tsv", sep = "", collapse = "") %>% system()
    tmp <- read.table("vf_b.tsv", header = FALSE, stringsAsFactors = FALSE,comment.char = "")
    colnames(tmp)<- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")

    results <- bind_rows(results,tmp)
  }
  results <- inner_join(results,annot, by = c("target" = "ID")) %>% separate(query,c("Genome","Protein"), sep = "#")
  return(results)
}


