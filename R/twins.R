#' Finds Twins into AccNET object
#'
#' Twins nodes are those nodes that share exactly the same edges. In the case
#' of an \emph{accnet} object twin are those genes/proteins that you can find
#' in the exactly the same set of genomes/samples.
#'
#' @param data An \emph{accnet} object
#'
#' @return A data.frame() with columns "Target", "Twin" and "Size".
#' \itemize{
#'   \item Target:  Column with gene/protein ID.
#'   \item Twin:  Twin label
#'   \item Size:  Twin group size (number of genes/proteins size)
#' }
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
twins <- function(data)
{
  if(!is(data,"accnet"))
  {
    stop("data must be an Accnet object")
  }

  data$matrix %>%
    column_to_rownames("Source") %>%
    t() %>%
    as.data.frame() %>%
    unite(pattern,1:nrow(data$matrix),sep = "") %>%
    rownames_to_column("Target") %>%
    group_by(pattern) %>%
    mutate(Freq = n(),Index = group_indices()) %>%
    ungroup() %>%
    mutate(Twin = ifelse(Freq>1,paste0("Twin",Index),"Single")) %>%
    select(Target,Twin,Size=Freq) %>%
      return()
}
