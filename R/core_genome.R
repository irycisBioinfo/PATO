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
#'
#' @return A core_genome object (a data.frame with two columns: fasta header and sequence)
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#' @import foreach
#' @import doParallel
#' @import parallel
#'
core_genome <- function(data, type, n_cores)
{
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

  if(sum(grep("avx2",system("cat /proc/cpuinfo",intern = TRUE),ignore.case = TRUE)))
  {
    mmseqPah = system.file("mmseqs.avx2", package = "pato")
  }else{
    mmseqPah = system.file("mmseqs.sse41", package = "pato")
  }

  blastParser = system.file("blast_parser.pl", package = "pato")

  if(file.exists("all.mmseq") & file.exists("all.cluster.index"))
  {

    data$table <- data$table %>% as_tibble()
    nGenomes  <-  data$table  %>% distinct(Genome_genome) %>% count()

    table <-  data$table %>%
      distinct() %>%
      group_by(Prot_prot) %>%
      mutate(Nprot = n_distinct(Genome_prot), Ngenomes = n_distinct(Genome_genome)) %>% ungroup() %>%
      filter(Ngenomes ==  nGenomes$n & Nprot == nGenomes$n)


    lookup <- data.table::fread("all.mmseq.lookup", sep = "\t") %>% as_tibble()
    colnames(lookup) = c("ID","head","value")

    core_table <-  table %>%
      unite(head, Prot_genome, Prot_prot, sep = "#", remove = FALSE) %>%
      select(head) %>%
      distinct() %>% inner_join(lookup)

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


     paste(mmseqPah," createsubdb IDS.txt all.cluster all.cluster.subset -v 0", sep = "", collapse = "") %>%
       system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
     paste(mmseqPah," createseqfiledb all.mmseq all.cluster.subset subset  -v 0", sep = "", collapse = "") %>%
       system(.,intern = F,ignore.stdout = T, ignore.stderr = T)
     paste(mmseqPah," result2flat all.mmseq all.mmseq subset subset.fasta", sep = "", collapse = "") %>%
       system(.,intern = F,ignore.stdout = T, ignore.stderr = T)

     paste("split -a 4 --numeric-suffixes=1 -l ",(nGenomes$n*2)+1," subset.fasta core_") %>%
       system(.,intern = F,ignore.stdout = T, ignore.stderr = T)


    system("mv core* f_core")
    system("sed -i '1d' ./f_core/*",intern = F,ignore.stdout = T, ignore.stderr = T)
    # seqs <-  data.frame()

    #for (i in(dir("./f_core")))

    #cl <- makeCluster(n_cores)
    #registerDoParallel(cl)

    #seqs <- foreach( i = dir("./f_core"), .combine = "rbind") %dopar%
    seqs = data.frame()
    for(i in dir("./f_core"))
    {
      #cat(i,'\n')
      system(paste("head -2 ./f_core/",i," > ./f_core/",i,".ref.fasta",sep = "",collapse = ""))
      system(paste("grep '>' ./f_core/",i," > ./f_core/headers_",i,sep = "",collapse = ""))
      l = system(paste("head -2 ./f_core/",i,"| tail -1 |wc -m ",sep = "",collapse = ""), intern =T)

      l = as.numeric(l)*10
      if(type =="nucl")
      {
        paste("blastn -query ./f_core/",i,".ref.fasta -subject ./f_core/",
              i,
              " -outfmt 4 -max_hsps 1 -out ./f_core/",i,".blast -line_length ",
              l,
              " -num_alignments ",
              nGenomes$n+10,
              sep = "", collapse = "") %>%
          system()
      }else if(type=="prot")
      {
        paste("blastp -query ./f_core/",i,".ref.fasta -subject ./f_core/",
              i,
              " -outfmt 4 -max_hsps 1 -out ./f_core/",i,".blast -line_length ",
              l,
              " -num_alignments ",
              nGenomes$n+10,
              sep = "", collapse = "") %>%
          system()
      }else{
        stop("type must be 'prot' or 'nucl'")
      }
      paste("perl ",blastParser," ./f_core/",i,".blast ./f_core/headers_",i," > ./f_core/",i,".aln", sep = "", collapse = "") %>%
        system()

      tmp <- readLines(paste("./f_core/",i,".aln",sep = "",collapse = ""))
      seqs <- bind_rows(seqs,data.frame(Head = tmp[seq(1,length(tmp)-1,2)],Seq = tmp[seq(2,length(tmp),2)]))

#
#       if(type=='nucl')
#       {
#         individual_aln[[i]] <- read.FASTA(paste("./f_core/",i,".aln",sep = "",collapse = ""), type = "DNA")
#       }else{
#         individual_aln[[i]] <- read.FASTA(paste("./f_core/",i,".aln",sep = "",collapse = ""), type = "AA")
#       }


    }
    #stopCluster(cl)

    print(colnames(seqs))
    seqs <- seqs %>%
      separate(Head,c("Genomes","Prot"),sep="#") %>%
      group_by(Genomes) %>%
      summarise(Seq = paste(Seq,sep = "",collapse = ""))

    print("Number of hard core-genes (100% presence split paralogous):")
    print(table %>% select(Prot_prot) %>% distinct() %>% nrow())

    results <- list(core_genome = seqs)
    class(results) <-  append(class(results),"core_genome" )
    return(results)

  } else{
    stop("You must execute mmseqs() before core_genome()")
  }






}
