#' HeatMap of Annotation (VF and/or AbR)
#'
#' Creates a heatmap with the annotation results. User can
#' filter the results by identity and/or evalue. This function
#' uses \emph{pheatmap} instead of \emph{heatmap} funtion if
#' \emph{pheatmap} is installed.
#'
#' @param data Data frame result from annotate.
#' @param min_identity Minimun identity to show in figure.
#' @param max_evalue Maximun Evalue to show in figure.
#'
#' @export
#'
#' @examples
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import dtplyr
#'
heatmap_of_annotation <- function(data, min_identity = 0.95, max_evalue = 1e-25)
{

  tmp <-  data %>%
    group_by(Genome,Protein) %>%
    summarise_all(first) %>%
    filter(pident >= min_identity) %>%
    filter(evalue <= max_evalue) %>%
    group_by(Genome,Gene) %>%
    summarise(pident = max(pident))%>%
    spread(Gene,pident, fill = 0) %>%
    column_to_rownames("Genome") %>% as.matrix()

    if(sum(grepl("pheatmap",installed.packages())))
    {
      pheatmap::pheatmap(tmp,clustering_method = "average",clustering_distance_rows = "binary", clustering_distance_cols = "binary")
    }else{
      heatmap(tmp)
    }
}
