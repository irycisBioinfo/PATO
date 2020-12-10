#' Title
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
maf_to_data.frame <- function(file)
{

  coords = system(paste("grep 'ref.' ",file," | cut -f 2,3,4,5,6 --delimiter=' '", sep = "",collapse = ""),intern = T)


  coords = coords %>%  as_tibble() %>%
    rename(all = 1) %>%
    separate(all,c("src","start","size","strand","srcSize"), sep = "\\s") %>%
    mutate(start = as.numeric(start)) %>%
    mutate(size = as.numeric(size)) %>%
    mutate(srcSize = as.numeric(srcSize)) %>%
    mutate(end = start+size)

  output = gsub(".maf",".bed",file)

  write.table(file = output,x = coords %>% select(src,start,end), col.names = F, row.names = F, quote = F, sep = "\t")

}
