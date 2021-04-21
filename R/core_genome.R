#' Core-Genome Alignment
#'
#' Find and creates a core-genome alignment. Unlike \emph{core_plots()}
#' this function find the hard core-genome (genes presence in 100% of genomes
#' and without repetitions). The function takes a \emph{mmseqs()} output, so the deffinition of
#' the paralogous genes of the core-genome (similarity, coverage and/or e-value)
#' depends on the \emph{mmseqs()} parameters.
#'
#' The function perform a pseudo-msa per each paralog using the function
#' \emph{result2msa} of \emph{mmseqs2}.This approach is much faster than
#' classical MSA (clutal, mafft or muscle) but is less accurate. Taking into account
#' that most of the phylogenetic inference software only takes variant columns with
#' no insertions or deletion, there are not to many difference in the final phylogenetic trees.
#'
#' \emph{core_genome()} can build a core-genome alingment of thusands of genomes in minutes.
#'
#' @param data An \emph{mmseqs} object
#' @param type Type of sequence 'nucl' or 'prot'
#' @param n_cores Number of computer core to use
#' @param methos \emph{fast (based on blast)} or \emph{accurate (based on mafft)}
#'
#' @return A core_genome object (a data.frame with two columns: fasta header and sequence)
#' @export
#'
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#' @import foreach
#' @import doParallel
#' @import parallel
#'
core_genome <- function(data, type, n_cores, method = "fast")
{
  if(Sys.which("perl")=="")
  {
    stop("This function needs perl to work. Please check that Perl is installed and in the PATH")
  }
  if(method == "fast" & Sys.which("blastn")=="")
  {
    stop("This function needs NCBI blastn to work. Please check that NCBI blast+ is installed and in the PATH")
  }
  if(method =="accurate" & Sys.which("mafft")=="")
  {
    stop("This function needs mafft to work. Please check that mafft is installed and in the PATH")
  }

  if(!is(data,"mmseq"))
  {
    stop("Error, 'data' must be a mmseq object")
  }

  if(missing(type))
  {
    stop("type parameter must be provided (nucl or prot)")
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


    if(grepl("AVX2",system("sysctl -a | grep 'AVX2'", intern = T)))
    {
      mmseqPath = system.file("mmseqs.macos.avx2", package = "pato")
    }else{
      mmseqPath = system.file("mmseqs.macos.sse41", package = "pato")
    }
  }else{
    stop("Error, OS not supported.")
  }
  blastParser = system.file("blast_parser.pl", package = "pato")

  msa2table = system.file("msa2table.pl", package = "pato")

  origin_path = getwd()
  setwd(data$path)
  on.exit(setwd(origin_path))

  table = data$table %>% as_tibble()

  nGenomes  <-  table  %>% distinct(Genome_genome) %>% count() %>% as_tibble()

  table <-  table %>%
    distinct() %>%
    group_by(Prot_prot) %>%
    mutate(Nprot = n_distinct(Genome_prot), Ngenomes = n_distinct(Genome_genome)) %>% ungroup() %>%
    filter(Ngenomes ==  nGenomes$n & Nprot == nGenomes$n) %>% as_tibble()


  lookup <- data.table::fread("all.mmseq.lookup", sep = "\t",stringsAsFactors = F,col.names = c("ID","head","value")) %>% as_tibble()


  core_table <-  table %>%
    unite(head, Prot_genome, Prot_prot, sep = "#", remove = FALSE) %>%
    select(head) %>%
    distinct() %>% inner_join(lookup)

  #core_table <- core_table %>% as_tibble()

  if(nrow(core_table) >0)
  {
    core_table %>%
      select(ID) %>% write.table("IDS.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }else{
    stop("No core_genome found")
  }

  if(dir.exists("f_core"))
  {
    system("rm -r f_core")
  }
  if(length(list.files(".",pattern = "core"))>0)
  {
    file.remove(list.files(".",pattern = "core"),recursive = TRUE)
  }

  dir.create("f_core")


  paste(mmseqPath," createsubdb IDS.txt all.cluster all.cluster.subset -v 0", sep = "", collapse = "") %>%
    system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
  paste(mmseqPath," createseqfiledb all.mmseq all.cluster.subset subset  -v 0", sep = "", collapse = "") %>%
    system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
  paste(mmseqPath," result2flat all.mmseq all.mmseq subset subset.fasta", sep = "", collapse = "") %>%
    system(.,intern = F,ignore.stdout = T, ignore.stderr = T)

  paste("split -a 4 --numeric-suffixes=1 -l ",(nGenomes$n*2)+1," subset.fasta core_") %>%
    system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
  system("mv core* f_core")
  system("sed -i '1d' ./f_core/*",intern = F,ignore.stdout = T, ignore.stderr = T)


  if(type =="nucl")
  {

    cl <- makeCluster(n_cores)
    registerDoParallel(cl)

    seqs <- foreach( i = dir("./f_core"), .combine = "rbind") %dopar%
    {
      if(method =="fast")
      {
      system(paste("head -2 ./f_core/",i," > ./f_core/",i,".ref.fasta",sep = "",collapse = ""))
      system(paste("grep '>' ./f_core/",i," > ./f_core/headers_",i,sep = "",collapse = ""))
      l = system(paste("head -2 ./f_core/",i,"| tail -1 |wc -m ",sep = "",collapse = ""), intern =T)
      l = as.numeric(l)*10

      paste("blastn -task blastn -query ./f_core/",i,".ref.fasta -subject ./f_core/",
             i,
             " -outfmt 4 -max_hsps 1 -out ./f_core/",i,".blast -line_length ",
             l,
             " -num_alignments ",
             nGenomes$n+10,
             sep = "", collapse = "") %>% system()
        paste("perl ",blastParser," ./f_core/",i,".blast ./f_core/headers_",i," > ./f_core/",i,".aln", sep = "", collapse = "") %>%
          system()
        tmp <- microseq::readFasta(paste0("./f_core/",i,".aln"))
        tmp <- tmp %>% mutate(Header = gsub("^",">",Header)) %>% mutate(Sequence = toupper(Sequence))

      }else if(method =="accurate")
      {

        options(warn = -1)
        system(paste0("mafft --quiet --thread 1 ./f_core/",i," > ./f_core/",i,".aln" ))
        options(warn = 0)
        tmp <- microseq::readFasta(paste0("./f_core/",i,".aln"))
        tmp <- tmp %>% mutate(Header = gsub("^",">",Header)) %>% mutate(Sequence = toupper(Sequence))
      }

    }

    stopCluster(cl)

    }else if(type=="prot")
    {
        paste(mmseqPath," result2msa all.mmseq all.mmseq all.cluster.subset all.core", sep = "", collapse = "") %>%
          system(.,intern = F,ignore.stdout = T, ignore.stderr = T)

        paste(mmseqPath," result2flat all.mmseq all.mmseq all.core all.core.fasta --use-fasta-header 0", sep = "", collapse = "") %>%
          system(.,intern = F,ignore.stdout = T, ignore.stderr = T)

        print(paste("perl ",msa2table," all.core.fasta > all.core.fasta.tab", sep = "", collapse = ""))
        paste("perl ",msa2table," all.core.fasta > all.core.fasta.tab", sep = "", collapse = "") %>%
          system()

        seqs <- data.table::fread("all.core.fasta.tab",
            sep = "\t",
            stringsAsFactors = F,
            header = F,
            colClasses = c("character", "character"),
            col.names = c("Head", "Seq")
          ) %>% as_tibble()

    }


  print(colnames(seqs))
  colnames(seqs) <- c("Head","Seq")
  seqs <- seqs %>% as_tibble() %>%
    separate(Head,c("Genomes","Prot"),sep="#") %>%
    group_by(Genomes) %>%
    summarise(Seq = paste(Seq,sep = "",collapse = ""))

  print("Number of hard core-genes (100% presence split paralogous):")
  print(table %>% select(Prot_prot) %>% distinct() %>% nrow())

  results <- list(core_genome = seqs,path = data$path)
  class(results) <-  append(class(results),"core_genome" )
  return(results)


}
