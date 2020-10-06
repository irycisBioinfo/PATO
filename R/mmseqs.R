#' MMSeqs2 Orthologous clustering.
#'
#' MMseqs2 (Many-against-Many sequence searching) is a software suite to
#' search and cluster huge protein and nucleotide sequence sets.
#' MMseqs2 is open source GPL-licensed software (https://mmseqs.com).
#' This function is a wrapper for the original \emph{mmseqs} executable.
#' The function clusters the genes or proteins in ortholog clusters and save the
#' representative member of each one and its annotation.
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
#' \item 3: target seq. length needs be at least x% of query length
#' \item 4: query seq. length needs
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
#' @export
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
#' @import foreach
#' @import doParallel
#'
mmseqs <- function(file_list, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores, cov_mode = 0, cluster_mode = 0)
{

  if(missing(n_cores))
  {
    n_cores = detectCores()-1
  }

  proc_cpu = readLines("/proc/cpuinfo")

  if(sum(grep("avx2",proc_cpu,ignore.case = TRUE)))
  {
    mmseqPath = system.file("mmseqs.avx2", package = "pato")
  }else{
    mmseqPath = system.file("mmseqs.sse41", package = "pato")
  }


  system("rm *.rnm")
  system("rm -r tmpDir")
  system("rm all*")
  system("rm commands.txt")

  for (i in file_list[,1])
  {
    if(grepl("gz",i[1]))
    {
      write(paste("zcat ",i," | perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' >> all.rnm \n", collapse = "",sep = ""),file = "commands.txt", append = T)
    }else{
      write(paste("perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' ",i," >> all.rnm \n", collapse = "",sep = ""),file = "commands.txt", append = T)

    }
  }

  system(paste(Sys.getenv("SHELL")," commands.txt",collapse = "",sep = ""))



  cmd1 <- paste(mmseqPath," createdb all.rnm all.mmseq",sep = "",collapse = "")
  print(cmd1)
  system(cmd1)

  cmd2 <- paste(mmseqPath," linclust all.mmseq all.cluster tmpDir --threads ",n_cores,
                " -e ",evalue,
                " --min-seq-id ",identity,
                " -c ",coverage,
                "--cov-mode", cov_mode,
                "--cluster-mode",cluster_mode,
                sep = "",collapse = "")
  print(cmd2)
  system(cmd2)

  cmd3 <- paste(mmseqPath," createtsv all.mmseq all.mmseq all.cluster all.cluster.tsv",sep = "",collapse = "")
  print(cmd3)
  system(cmd3)

  cmd4 <- paste(mmseqPath," result2repseq all.mmseq all.cluster all.representatives",sep = "",collapse = "")
  print(cmd4)
  system(cmd4)

  cmd5 <- paste(mmseqPath," result2flat all.mmseq all.mmseq all.representatives all.representatives.fasta --use-fasta-header",sep = "",collapse = "")
  print(cmd5)
  system(cmd5)

  system("grep '>' all.representatives.fasta | sed 's/>//' | sed 's/ /~~/' | sed 's/\t/ /g' > AnnotFile.tsv");

  Annotation <- data.table::fread("AnnotFile.tsv", header = FALSE, sep = "\t") %>% as_tibble() %>%
  separate(V1, c("Genome","kk"), sep = "#") %>%
    separate("kk",c("Prot_prot","Annot"), sep= "~~") %>%
    select(Prot_prot,Annot)
  mmseqs.raw <- fread("all.cluster.tsv", header = FALSE, sep = "\t") %>% as_tibble()
  colnames(mmseqs.raw) = c("Prot","Genome")
  mmseqs.raw <- mmseqs.raw %>%
    separate(Prot,c("Prot_genome","Prot_prot"),sep = "#") %>%
    separate(Genome,c("Genome_genome","Genome_prot"), sep = "#")

  results <- list(table = mmseqs.raw, annot = Annotation)
  class(results) <- "mmseq"
  return(results)

}

