#' Creates a matrix of sharedness
#'
#' Count the number of genes/proteins that share the samples among them.
#'
#' @param data An \emph{accnet} object.
#'
#' @return A matrix with the count of the genes/proteins shared between
#' each pair of samples.
#' @export
#'
#' @examples
sharedness <- function(data)
{
  if(is(data,"accnet"))
  {
    if(is(data,"pangenome"))
    {

      tmp <- data$list %>% mutate(B = as.numeric(Weight >0)) %>%
        group_by(Source,Target) %>%
        summarise(B = max(B)) %>%
        spread(Target,B, fill = 0) %>%
        column_to_rownames("Source") %>%
        as.matrix()

    }else{
      tmp <- data$matrix %>% column_to_rownames("Source") %>% as.matrix()
    }



  }else if(is(data,"mmseq"))
  {
    tmp <- data$table %>% as_tibble() %>% ungroup() %>%
      select(Prot_prot, Genome_genome) %>%
      distinct() %>%
      mutate(presence =1) %>%
      spread(Prot_prot, presence, fill = 0) %>%
      column_to_rownames("Genome_genome") %>%
      as.matrix()
  }else{
    stop("Error 'data' must be accnet object")
  }

  result <- crossprod(t(tmp))

  #rownames(result) <- rownames(tmp)
  #colnames(result) <- rownames(tmp)
  return(result)

}
