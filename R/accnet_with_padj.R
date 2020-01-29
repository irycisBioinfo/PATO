#' Creates a AccNET network using the enrichment analysis
#'
#' Creates and improved network using the adjusted p-value as edge-weight.
#' The network could be export to gephi (\code{export_to_gephi}) for a correct
#' visualitation
#'
#' @param data \code{accnet_enrichment_analysis} result.
#'
#' @return Return and \code{accnet} object with the adjusted p-value as edge-weigth
#' @export
#'
#' @seealso \code{\link{accnet_enrichment_analysis}}, \code{\link{export_to_gephi}}
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table
#'
accnet_with_padj <- function(data)
{
	if(!is(data,"accnet_enr"))
  {
    stop("Data must be an 'accnet_enr' object")
  }
  list <- data %>%  select(Source,Target,Weight = padj) %>% mutate(Weight = 2-Weight)
  matrix <- data %>%  select(Source,Target,Weight = padj) %>% mutate(Weight = 2-Weight) %>% spread(Target,Weight,fill = 0)
  Annot <- data %>% select(ID = Target,Annot) %>% distinct()

  results <- list(list = list, matrix = matrix, annot = Annot)
  class(results) <- "accnet"
  return(results)
}
