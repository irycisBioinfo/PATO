#' Classifier
#'
#' Classifier take a list of genome files (nucleotide o protein) and identify
#' the most similar specie to each file.
#' Classsifier uses all reference and representative genomes from NCBI Refseq
#' database and search, using mash, the best hit for each genome file in the input list.
#'
#' @param file_list Data frame with the full path to the genome files (gene or protein multi-fasta) or \emph{gff_list} object.
#' @param n_cores Number of cores to use.
#' @param type Type of sequence 'nucl' (nucleotides), 'prot' (aminoacids) or 'wgs' for whole genome sequence (only with gff_list objects)
#' @param max_dist Maximun distance to report (1-Average Nucleotide Identity). Usually all species have 0.05 distance among all each memebers
#'
#' @return Classifier returns a data.frame with the best hit for each input genome.
#' @export
#'
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @importFrom  data.table fread

classifier <- function(file_list, n_cores, type ="nucl", max_dist = 0.06)
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
    }else if(type =="wgs"){
      file_list = dir(paste(file_list$path,"/fna/",sep = "", collapse = ""),full.names = T) %>% as_tibble()
    }
  }

  if(missing(n_cores))
  {
    n_cores <- detectCores()-1
  }

  if(type == "prot")
  {
    reference <- system.file("DB/ReferenceProt.msh",package = "pato")
  }else if(type == "nucl" | type == "wgs")
  {
    reference <- system.file("DB/ReferenceNucl.msh",package = "pato")
  }else{
    stop("Error in type options. Only prot or nucl options are allowed")
  }

  if(grepl('linux',Sys.getenv("R_PLATFORM"))) ## Linux
  {
    mashPath = system.file("mash",package = "pato")
  }else if(grepl('apple',Sys.getenv("R_PLATFORM"))){ ##MacOS
    mashPath = system.file("mash.macos",package = "pato")
  }else{
    stop("Error, OS not supported.")
  }


  folderName = paste(getwd(),"/",md5(paste(file_list[,1], sep = "",collapse = "")),"_mash",sep = "",collapse = "")

  if(!file.exists(paste(folderName,"/all.msh", sep = "",collapse = "")))
  {
    dir.create(folderName)
    write.table(file_list[,1],paste(folderName,"/input_mash.txt",sep = "",collapse = ""),
                quote = F, col.names = FALSE, row.names = FALSE)

    if(type == "prot")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -l ",folderName,"/input_mash.txt"," -a -o ",folderName,"/all.msh", sep = "", collapse = "")
    }else if(type =="nucl" | type =="wgs")
    {
      cmd1 <- paste(mashPath," sketch -p ",n_cores," -l ",folderName,"/input_mash.txt"," -o ",folderName,"/all.msh", sep = "", collapse = "")
    } else{
      stop("Error in type options. Only prot or nucl options are allowed")
    }

    system(cmd1)
  }

  cmd3 <- paste(mashPath," dist -p ",n_cores," -d ",max_dist," ",reference," ",folderName,"/all.msh > ",folderName,"/classification.tab", sep = "", collapse = "")
  system(cmd3)


  class.table <- data.table::fread(paste(folderName,"/classification.tab",sep = "",collapse = ""),
                                         header = F) %>% as_tibble()
  colnames(class.table) <- c("Source","Target","Dist","pvalue","sketch")


  header <- data.table::fread(system.file("DB/header_class.tsv",package = "pato"), sep = "\t", header = T, stringsAsFactors = F)



  class.table %>%
    mutate(Target = basename(Target)) %>%
    separate(Source,c("kk","acc1","acc2"), sep = "_") %>%
    unite(assembly_accession,kk,acc1,sep = "_") %>%
    mutate(Dist = 1-Dist) %>%
    left_join(header) %>%
    select(Target,Similarity = Dist,organism_name,species_taxid,assembly_accession,refseq_category) %>%
    group_by(Target) %>% slice_max(order_by = Similarity,n = 1) %>% ungroup() %>%
    return()


}
