#' Title
#'
#' @param mmseq
#' @param accnet
#' @param min_size
#'
#' @return
#' @export
#'
#' @examples
dn_ds <- function(mmseq,accnet,min_size,n_cores,mode)
{
  if(Sys.which("perl")=="")
  {
    stop("This function needs perl to work. Please check that Perl is installed and in the PATH")
  }
  if(Sys.which("blastn")=="")
  {
    stop("This function needs NCBI blastn to work. Please check that NCBI blast+ is installed and in the PATH")
  }

  if(!is(mmseq,"mmseq"))
  {
    stop("Error, 'mmseq' must be a mmseq object")
  }
  if(!is(accnet,"accnet"))
  {
    stop("Error, 'accnet' must be a accnet object")
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
  blastParser = system.file("blast_parser_dnds.pl", package = "pato")
  cluster_parser = system.file("cluster_parser.pl", package ="pato")


  if(dir.exists(paste0(mmseq$path,"/folder_0")))
  {
    system(paste0("rm -r ",mmseq$path,"/folder_*"))
  }


  origin_path = getwd()
  setwd(mmseq$path)
  on.exit(setwd(origin_path))
  on.exit(stopCluster(cl), add = TRUE)


  system(paste0(mmseqPath, " createseqfiledb all.mmseq all.cluster dnds"))
  system(paste0(mmseqPath," result2flat all.mmseq all.mmseq dnds dnds.fasta"))
  system(paste0("perl ",cluster_parser))




  folders <- dir(pattern = "folder_",recursive = F)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  results <- data.frame()
  for (f in folders)
  {
    setwd(f)
    system("grep -c '>' family_* | sed 's/:/\t/' > sizes.tsv")
    sizes <- read.table("sizes.tsv", sep = "\t", header = F)
    sizes <- sizes %>% rename(Family = 1, Size = 2) %>% filter(Size > min_size)

    if(mode =="fast"){
      dnds_ratios <- foreach( i = 1:nrow(sizes), .combine = "rbind") %dopar%
        {

          name = system(paste0("head -1 ",f,"/",sizes$Family[i]), intern = T)
          system(paste0("sed -i '$ d' ",f,"/",sizes$Family[i]))
          system(paste("head -2 ",f,"/",sizes$Family[i]," > ",f,"/",sizes$Family[i],".ref.fasta",sep = "",collapse = ""))
          l = system(paste("head -2 ",f,"/",sizes$Family[i],"| tail -1 |wc -m ",sep = "",collapse = ""), intern =T)
          system(paste("grep '>' ",f,"/",sizes$Family[i]," > ",f,"/headers_",sizes$Family[i],sep = "",collapse = ""))
          l = as.numeric(l)*10
          paste("blastn -task blastn -query ",f,"/",sizes$Family[i],".ref.fasta -subject ",f,"/",sizes$Family[i]," -outfmt 4 -max_hsps 1 -line_length ",
                l+10,
                " -num_alignments ",
                sizes$Size[i]," > ",f,"/",sizes$Family[i],".blast",
                sep = "", collapse = "") %>% system()
          paste("perl ",blastParser," ",f,"/",sizes$Family[i],".blast ",f,"/headers_",sizes$Family[i]," > ",f,"/",sizes$Family[i],".aln", sep = "", collapse = "") %>%
            system()


          F = read.dna(paste0(f,"/",sizes$Family[i],".aln"),format = "fasta")
          F = unique.matrix(F)

          if(nrow(F)>=2)
          {
            A = dnds(F,codonstart = 1) %>% as.matrix()
            diag(A) = NaN
            A[is.infinite(A)] = max(A[!is.infinite(A)])
            A[A<0] = NaN
            M = mean(A[!is.infinite(A)],na.rm = T)
            data.frame(Target = name, av_dnds = M)


          }else{
            data.frame(Target = name, av_dnds = NaN)
          }

        }

    }else if(mode =="accurate")
    {
      dnds_ratios <- foreach( i = 1:nrow(sizes), .combine = "rbind") %dopar%
        {
          name = system(paste0("head -1 ",f,"/",sizes$Family[i]), intern = T)
          system(paste0("sed -i '$ d' ",f,"/",sizes$Family[i]))
          options(warn = -1)
          system(paste0("mafft --quiet --thread 1 ",f,"/",sizes$Family[i]," > ",f,"/",sizes$Family[i],".aln"))
          options(warn = 0)

          F = read.dna(paste0(f,"/",sizes$Family[i],".aln"),format = "fasta")
          F = unique.matrix(F)

          if(nrow(F)>=2)
          {
            A = dnds(F,codonstart = 1) %>% as.matrix()
            diag(A) = NaN
            A[is.infinite(A)] = max(A[!is.infinite(A)])
            A[A<0] = NaN
            M = mean(A[!is.infinite(A)],na.rm = T)
            data.frame(Target = name, av_dnds = M)


          }else{
            data.frame(Target = name, av_dnds = NaN)
          }
        }
    }
    results <- bind_rows(results,dnds_ratios)
    setwd(mmseq$path)

  }

  accnet$annot <- results  %>%
    separate(Target,c("Genome","ID"), sep = "#") %>%
    select(ID,av_dnds) %>%
    right_join(accnet$annot)

  return(accnet)

}
