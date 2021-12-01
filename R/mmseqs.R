#' MMSeqs2 Orthologous clustering.
#'
#' MMseqs2 (Many-against-Many sequence searching) is a software suite to
#' search and cluster huge protein and nucleotide sequence sets.
#' MMseqs2 is open source GPL-licensed software (https://mmseqs.com).
#' This function is a wrapper for the original \emph{mmseqs} executable.
#' The function clusters the genes or proteins in ortholog clusters and save the
#' representative member of each one and its annotation.
#'
#' @param file_list Data frame with the full path to the genome files (gene or protein multi-fasta) or a \emph{gff_list} object.
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
#' @param type 'only with \emph{gff_list} objects. Specify what kind of data you want to use: 'nucl' or 'prot'
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
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#' @import foreach
#' @import doParallel
#' @import openssl
#'
mmseqs <- function(file_list, coverage = 0.8, identity = 0.8, evalue = 1e-6, n_cores, cov_mode = 0, cluster_mode = 0, type)
{

  if(is(file_list,"gff_list"))
  {
       if(missing(type))
       {
         stop("type must be declared for gff_list objects")
       }else if(type == "prot")
       {
         file_list = dir(paste(file_list$path,"/faa/",sep = "", collapse = ""),full.names = T) %>% as_tibble()
       }else if(type == "nucl")
       {
         file_list = dir(paste(file_list$path,"/ffn/",sep = "", collapse = ""),full.names = T) %>% as_tibble()
       }else{
         stop("Error in 'type' parameter. Option not recognized")
       }
  }else{
    file_list <- as.data.frame(file_list)
  }

  if(missing(n_cores))
  {
    n_cores = detectCores()-1
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

      mmseqPath = system.file("mmseqs.macos", package = "pato")
  }else{
    stop("Error, OS not supported.")
  }

  folderName = paste(getwd(),"/",md5(paste(file_list[,1], sep = "",collapse = "")),"_mmseq",sep = "",collapse = "")

  if(!dir.exists(folderName))
  {
    dir.create(folderName,)
  }

  # system("rm *.rnm")
  # system("rm -r tmpDir")
  # system("rm all*")


  if(!file.exists(paste(folderName,"/all.rnm",sep = "",collapse = "")))
  {
    if(file.exists("commands.txt"))
    {
      file.remove("commands.txt")
    }

    for (i in file_list[,1])
    {
      if(grepl("gz",i[1]))
      {
        write(paste("zcat ",i," | perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' >> ",folderName,"/all.rnm \n", collapse = "",sep = ""),
              file = "commands.txt",
              append = T)
      }else{
        write(paste("perl -pe 's/>/$&.\"",basename(i),"\".\"#\".++$n.\"|\"/e' ",i," >> ",folderName,"/all.rnm \n", collapse = "",sep = ""),
              file = "commands.txt",
              append = T)
      }
    }
    system(paste(Sys.getenv("SHELL")," commands.txt",collapse = "",sep = ""))
  }



  origin_path = getwd()
  setwd(folderName)
  on.exit(setwd(origin_path))



  if(!file.exists("all.mmseq"))
  {
    cmd1 <- paste(mmseqPath," createdb all.rnm all.mmseq",sep = "",collapse = "")
    print(cmd1)
    system(cmd1)
  }

  if(!file.exists("all.cluster"))
  {
    cmd2 <- paste(mmseqPath," linclust all.mmseq all.cluster . --threads ",n_cores,
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

    if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
    {
      system("sed -i 's/#/\t/g' all.cluster.tsv")
      mmseqs.raw <- fread("all.cluster.tsv", header = FALSE, sep = "\t")
    }else{
      system("perl -pe 's/#/\t/g' all.cluster.tsv > tmp")
      mmseqs.raw <- fread("tmp", header = FALSE, sep = "\t")
    }

    colnames(mmseqs.raw) = c("Prot_genome","Prot_prot","Genome_genome","Genome_prot")


  }


  results <- list(table = mmseqs.raw, annot = Annotation, path = folderName)
  class(results) <- append(class(results),"mmseq")
  return(results)

}

