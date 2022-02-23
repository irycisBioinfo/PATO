#' Screening
#'
#' @param data A \code{mmseq} object
#' @param type user must be specified if the data set is nucleotide, wgs or protein.
#' @param database A vector with the query databases:
#' \itemize{
#' \item \emph{AbR}: Antibiotic resistance database (ResFinder)
#' \item \emph{VF_A} VFDB core dataset (genes associated with experimentally verified VFs only)
#' \item \emph{VF_B} VFDB full dataset (all genes related to known and predicted VFs)
#' \item \emph{bacmet} BacMet dataset (genes associated with Biocide and Metal resistance)
#' }
#' @param query "all" or "accessory". It perform the DB from whole
#' protein dataset or just from the accessory

#'
#' @return
#' @export
#'
#' @examples
screening <- function(data, type = "nucl", database =c("AbR","VF_A","VF_B"), query = "all")
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

  if(!is(data,"mmseq"))
  {
    stop("Error data must be a mmseq object")
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

  n_cores = detectCores()

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
  }else{
    stop("Error in data type selection: please specify 'nucl' or 'prot'")
  }

  origin <- getwd()
  setwd(data$path)
  on.exit(setwd(origin))






}
