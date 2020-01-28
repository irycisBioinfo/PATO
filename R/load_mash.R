#' Load MASH
#'
#' Load an output from MASH. The file could be the matrix distance
#' (-t option in mash command line) or the standart table.
#'
#' @param file Path to mash output file
#' @param format \emph{list} or \emph{matrix} type.
#'
#' @return A \emph{mash} object
#'
#' @export
#'
#' @examples
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#' @import data.table


load_mash <- function(file, format = "matrix")
{

  if(format == "matrix")
  {
    mash.matrix <- fread(file, header = TRUE)
    colnames(mash.matrix) <- gsub("#","",colnames(mash.matrix))

    mash.list <- mash.matrix %>%
      mutate(Genome = basename(query)) %>%
      select(-query)%>%
      gather(Target,Dist, -Genome) %>%
      rename(Source = Genome) %>%
      mutate(Target = basename(Target)) %>% as.data.table()

    mash.matrix <- mash.matrix %>%
      mutate(Genome = basename(query)) %>%
      select(-query) %>%
      column_to_rownames("Genome") %>%
      as.matrix()

    colnames(mash.matrix) <- rownames(mash.matrix)
    results <- list(matrix = mash.matrix, table = mash.list)
    class(results) <- "mash"
    return(results)
  }else if(format =="list")
  {
    mash.list <- data.table::fread(file,header = FALSE)
    colnames(mash.list) <- c("Source","Target","Dist","pvalue","sketch")
    mash.list <- mash.list %>% select(-pvalue, -sketch)
    mash.matrix <- mash.matrix %>%
      spread(Target,Dist, fill = 1) %>%
      column_to_rownames("Target") %>% as.data.table()
    results <- list(matrix = mash.matrix %>% as.matrix(), table = mash.list %>% as_tibble())
    class(results) <- "mash"
    return(results)
  }else{
    stop("Error in format argument. Please specify 'matrix' or 'list' format")
  }
}
