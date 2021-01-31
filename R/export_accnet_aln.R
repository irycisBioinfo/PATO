#' Export accesory genome binary multi-alignment.
#'
#' This function write the presence/abscence matrix into a binary alignment file
#' that can be use with RaxML or iqtree to build a phylogenetic tree.
#'
#' Some phylogenetic tree inference tool allow as input binary aligments. We recommend
#' to use RAxML-NG or iQtree to inference an accessory genome tree.
#'
#' @param accnet An \emph{accnet} object
#' @param file Output filename
#' @param min_freq Minimun frequency of a gene/proteins to include in the alignment
#'
#' @export
#'
#' 
#'
export_accnet_aln <- function(accnet,file, min_freq =3)
{
  if(file.exists(file))
  {
    file.remove(file)
  }
  if(is(accnet,"accnet"))
  {
    acc_bin <- accnet$list %>%
      filter(degree >= min_freq) %>%
      mutate(value ="A") %>%
      select(Source,Target,value) %>%
      distinct() %>%
      spread(Target,value, fill = "C") %>%
      unite(Seq, -Source, sep = "")


    acc_bin$Source  <-  gsub("^",">",acc_bin$Source)

    for(i in 1:nrow(acc_bin))
    {
      write(acc_bin$Source[i],file,append = TRUE)
      write(acc_bin$Seq[i],file, append = TRUE)
    }
  }
}
