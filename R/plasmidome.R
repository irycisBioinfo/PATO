#' Plasmidome: Extract the plasmid sequences from a set of genomes
#'
#' This function uses \emph{mlplasmids} to find the contigs belong to plasmids.
#' The package \emph{mlplasmids} must be installed in your computer (see datails).
#' Nowadays \emph{mlplasmdis} only process Enterococcus faecium, Escherichia coli
#' or Klebsiella pneumoniae genomes
#'
#' The input genomes must be in \emph{gff_list} format (see \emph{load_gff_list()} function)
#' Plasmidome find the contigs that belong to plasmids and creates a new GFF structure
#' with all the features information (FAA, FNA, FFN and GFF).
#'
#' By default PATO does not install \emph{mlplasmids} because it is deposits in
#' its own gitlab repository (https://gitlab.com/sirarredondo/mlplasmids).
#'
#' The authors define \emph{mlplasmids} as:
#'
#' \cite{mlplasmids consists of binary classifiers to predict contigs either as
#' plasmid-derived or chromosome-derived. It currently classifies short-read
#' contigs as chromosomes and plasmids for Enterococcus faecium, Escherichia coli
#' and Klebsiella pneumoniae. Further species will be added in the future. The
#' provided classifiers are Support Vector Machines which were previously
#' optimized and selected versus other existing machine-learning techniques
#' (Logistic regression, Random Forest..., see also Reference). The classifiers
#' use pentamer frequencies (n = 1,024) to infer whether a particular contig
#' is plasmid- or chromosome- derived.}
#'
#' To install \emph{mlplasmids} you must have installed \emph{devtools} (usually
#' you have because PATO requires that you have devtools installed too). Then type:
#'
#' \code{devtools::install_git("https://gitlab.com/sirarredondo/mlplasmids")}
#'
#' @param gff_list A gff_list object (see \code{load_gff_list()})
#' @param specie mlplasmid accept as species Enterococcus faecium, Escherichia coli
#' or Klebsiella pneumoniae
#' @return
#' @export
#'
#'
#' @references Arredondo-Alonso et al., Microbial Genomics 2018;4 DOI 10.1099/mgen.0.000224
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom data.table fread
#' @import openssl
#' @import microseq
#'
#' @value return a \emph{gff_list} object


#'
plasmidome <- function(gff_list, specie)
{
  if(!nzchar(system.file(package = "mlplasmids")))
  {
    message("
    Error: mlplasmids package not found.

    This functior require mlplasmids package
    The mlplamids package is available through a gitlab repository

    Please use

    > install.packages(\"devtools\")
    > devtools::install_git(\"sirarredondo/mlplasmids\")

    to install the package.")
    stop()
  }else{
    require(mlplasmids)
  }

  if(missing(specie))
  {
    stop("Error: you must specified an specie: Enterococcus faecium, Escherichia coli o Klebsiella pneumoniae")
  }


  folderName <- gsub("_gffList","_plasmidome",gff_list$path)
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


  fastas <- dir(paste0(gff_list$path,"/fna/"),pattern = ".fna", full.names = T)
  gffs <- dir(paste0(gff_list$path,"/gff/"),pattern = ".gff", full.names = T)



  for(i in 1:length(fastas))
  {
    mlp <- plasmid_classification(fastas[i], species = specie)
    if(nrow(mlp)>0)
    {

      gff <- readGFF(gffs[i])
      fasta <- readFasta(fastas[i])

      fasta <- fasta %>% separate(Header, c("Seqid"), sep=" ", remove =F)
      mlp <- mlp %>% separate(Contig_name,c("Seqid"), sep=" ", remove =F)

      gff <- semi_join(gff,mlp)
      fasta <- semi_join(fasta,mlp) %>% select(-Seqid)
      pathName<- gsub(".fna","",basename(as.character(fastas[i])))

      writeGFF(gff,paste(folderName,"/gff/",pathName,".gff",sep = "",collapse = "") )                         ##Write the gff file
      writeFasta(fasta,paste(folderName,"/fna/",pathName,".fna",sep = "",collapse = ""))                      ## Write the fna file
      ffn_faa <- gff2fasta(gff %>% filter(Type == "CDS"),fasta)

      writeFasta(ffn_faa, paste(folderName,"/ffn/",pathName,".ffn",sep = "",collapse = "")) ## Write the ffn file
      ffn_faa <- ffn_faa %>% mutate(Sequence = microseq::translate(Sequence))
      writeFasta(ffn_faa , paste(folderName,"/faa/",pathName,".faa",sep = "",collapse = "")) ## Write the faa file
    }
  }

  results <-  list(path = folderName, files = gff_list$input_files)
  class(results) <- append(class(results),"gff_list")
  saveRDS(results,file = paste(folderName,"/gffObject.Rdata",sep = "",collapse = ""))
  return(results)
}
