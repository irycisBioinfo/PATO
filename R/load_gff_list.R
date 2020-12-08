#' Load GFF3 file list
#'
#' This function load and process a list of GFF3 files. The function extract CDS in nucleotide and
#' in aminoacids. It creates four folder with the GFF anotation, the FNA whole-genome sequence,
#' the FFN for fasta features nucleotide and the FAA with the proteome function.
#'
#' Load GFF3 files is slower than load multiple-fasta files (nucleotide or aminoacids) so for large
#' datasets we recommend to load fasta files instead of GFF3.
#'
#' @param file
#'
#' @return
#' The function returns a \code{gff_list} object the can be used as input for other functions (mmseqs, mash)
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#' @import Biostrings
#' @import openssl
#' @import microseq
#' @import foreach
#' @import doParallel
#'
load_gff_list <- function(input_files, n_cores)
{

  if(missing(n_cores)){

    n_cores = detectCores()-1
  }

  cl <- makeCluster(n_cores)
  registerDoParallel(cl)

  input_files <- as_tibble(input_files) %>% rename(File = 1)
  folderName <- paste(getwd(),"/",md5(paste(input_files$File, sep = "",collapse = "")),"_gffList",sep = "",collapse = "")

  if(file.exists(paste(folderName,"/gffObject.Rdata", sep = "", collapse = "")))
  {
    return(readRDS(paste(folderName,"/gffObject.Rdata", sep = "", collapse = "")))
  }

  if(!dir.exists(folderName))
  {
    dir.create(folderName)
    dir.create(paste(folderName,"/gff",sep = "",collapse = ""))
    dir.create(paste(folderName,"/fna",sep = "",collapse = ""))
    dir.create(paste(folderName,"/faa",sep = "",collapse = ""))
    dir.create(paste(folderName,"/ffn",sep = "",collapse = ""))
  }

  foreach (i = 1:nrow(input_files)) %dopar%
  {

    tempFile = tempfile(pattern = "gffTMP")
    print(as.character(input_files$File[i]))
    system(paste("csplit ",as.character(input_files$File[i])," -f ",tempFile," /#FASTA/"),ignore.stdout = T)

    pathName<- gsub(".gff","",basename(as.character(input_files$File[i])))

    gff <- readGFF(paste(tempFile,"00",sep = "",collapse = ""))
    writeGFF(gff,paste(folderName,"/gff/",pathName,".gff",sep = "",collapse = "") )                         ##Write the gff file

    fasta <- readFasta(paste(tempFile,"01",sep = "",collapse = ""))
    writeFasta(fasta,paste(folderName,"/fna/",pathName,".fna",sep = "",collapse = ""))                      ## Write the fna file


    ffn_faa <- gff2fasta(gff %>% filter(Type == "CDS"),fasta)

    writeFasta(ffn_faa, paste(folderName,"/ffn/",pathName,".ffn",sep = "",collapse = "")) ## Write the ffn file
    ffn_faa <- ffn_faa %>% mutate(Sequence = translate(Sequence))
    writeFasta(ffn_faa , paste(folderName,"/faa/",pathName,".faa",sep = "",collapse = "")) ## Write the faa file
  }


  on.exit(file.remove(list.files(pattern = "gffTMP")))
  on.exit(stopCluster(cl), add = T)

  results <-  list(path = folderName, files = input_files)
  class(results) <- append(class(results),"gff_list")
  saveRDS(results,file = paste(folderName,"/gffObject.Rdata",sep = "",collapse = ""))
  return(results)

}
