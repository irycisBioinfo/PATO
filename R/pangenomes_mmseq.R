#' Final step of pangenomes homologous cluster
#'
#' This is an internal function for pagenomes workflow
#'
#' @param file_list Data frame with the full path to the genome files (gene or protein multi-fasta).
#' @param coverage Minimun coverage (length) to cluster.
#' @param identity Minimun Identity.
#' @param evalue Maximun Evalue.
#' @param n_cores Number of cores to use.
#' @param cluster_mode 	Cluster mode:\itemize{
#' \item{0: Setcover}
#' \item{1: connected component}
#' \item{2: Greedy clustering by sequence length}
#' \item{3: Greedy clustering by sequence length (low mem)}
#' }
#' @param cov_mode Coverage mode:\itemize{
#' \item 0: Coverage of query and target
#' \item 1: Coverage of target
#' \item 2: coverage of query
#' \item 3: target seq.length needs be at least x% of query length
#' \item 4: query seq.length needs
#' }
#'
#' @return Return a \emph{mmseq} object.
#'
#' @note A \emph{mmseq} object is a list of two elements.
#' First contains a data.table/data.frame with four columns (Prot_genome, Prot_Prot,
#' Genome_genome and Genome_Prot). This is the output of MMSeqs2 and described the clustering
#' of the input genes/proteins. First column referes to the genome that contain the
#' representative gene/protein of the cluster. Second one, is the representative protein of the
#' cluster (i.e. the cluster name). Third colum is the genome that contains the gene/protein of
#' the fourth column.
#'
#' In the second element we can find a data.frame/data.table with the original annotation of all
#' representative gene/protein of each cluster in two columns. The first one Prot_prot is
#' the same that the second one of the first element.
#'
#'
#' @references Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi: 10.1038/nbt.3988 (2017).
#' @references Steinegger M and Soeding J. Clustering huge protein sequence sets in linear time. Nature Communications, doi: 10.1038/s41467-018-04964-5 (2018).
#' @references Mirdita M, Steinegger M and Soeding J. MMseqs2 desktop and local web server app for fast, interactive sequence searches. Bioinformatics, doi: 10.1093/bioinformatics/bty1057 (2019)
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
pangenomes_mmseqs <- function(file_list, coverage, identity, evalue, n_cores, cov_mode, cluster_mode,folder)
{

  if(sum(grep("avx2",system("cat /proc/cpuinfo",intern = TRUE),ignore.case = TRUE)))
  {
    mmseqPah = system.file("mmseqs.avx2", package = "pato")
  }else{
    mmseqPah = system.file("mmseqs.sse41", package = "pato")
  }


  system("rm pang.all.* ")
  system("rm -r pangtmpDir")

  for (i in file_list[,1])
  {
    fields = strsplit(i,split = "\\/", perl = TRUE)

    cmd <- paste("sed 's/>/>",basename(i),"#/' ",folder,"/",i," >> ",folder,"/pang.all.rnm", collapse = "",sep = "")
    print(cmd)
    system(cmd)
  }

  print(cmd)


  cmd <- paste(mmseqPah," createdb ",folder,"/pang.all.rnm ",folder,"/pang.all.mmseq -v 0",sep = "",collapse = "")
  print(cmd)
  system(cmd)



  cmd <- paste(mmseqPah," linclust ",folder,"/pang.all.mmseq ",folder,"/pang.all.cluster ",folder,"/pangtmpDir  -v 0 --threads ",n_cores,
               " -e ",evalue,
               " --min-seq-id ",identity,
               " -c ",coverage,
               " --cov-mode ", cov_mode,
               " --cluster-mode ",cluster_mode,
               sep = "",collapse = "")
  print(cmd)
  system(cmd)


  cmd <- paste(mmseqPah," createtsv ",folder,"/pang.all.mmseq ",folder,"/pang.all.mmseq ",folder,"/pang.all.cluster ",folder,"/pang.all.cluster.tsv  -v 0",sep = "",collapse = "")
  print(cmd)
  system(cmd)



  cmd <- paste(mmseqPah," result2repseq ",folder,"/pang.all.mmseq ",folder,"/pang.all.cluster ",folder,"/pang.all.representatives  -v 0",sep = "",collapse = "")
  print(cmd)
  system(cmd)



  cmd <- paste(mmseqPah," result2flat ",folder,"/pang.all.mmseq ",folder,"/pang.all.mmseq ",folder,"/pang.all.representatives ",folder,"/pang.all.representatives.fasta  -v 0 --use-fasta-header",sep = "",collapse = "")
  print(cmd)
  system(cmd)


  system(paste("grep '>' ",
               folder,
               "/pang.all.representatives.fasta | sed 's/>//' | sed 's/ /~~/' | sed 's/\t/ /g' > ",
               folder,
               "/pang.AnnotFile.tsv", sep = "",collapse = ""))

  Annotation <- data.table::fread(paste(folder,"/pang.AnnotFile.tsv",sep="",collapse = ""), header = FALSE, sep = "\t")
  Annotation <- Annotation %>% separate(V1, c("Genome","kk"), sep = "#") %>%
    separate("kk",c("Prot_prot","Annot"), sep= "~~") %>%
    select(Prot_prot,Annot)
  mmseqs.raw <- fread(paste(folder,"/pang.all.cluster.tsv",sep="",collapse = ""), header = FALSE, sep = "\t")
  colnames(mmseqs.raw) = c("Prot","Genome")
  mmseqs.raw <- mmseqs.raw %>%
    separate(Prot,c("Prot_genome","Prot_prot"),sep = "#") %>%
    separate(Genome,c("Genome_genome","Genome_prot"), sep = "#")

  results <- list(table = mmseqs.raw, annot = Annotation)
  class(results) <- "mmseq"
  return(results)

}

