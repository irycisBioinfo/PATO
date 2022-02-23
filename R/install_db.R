#' Install PATO reference data bases
#'
#' This function download and install the reference data bases for the functions
#' \emph{screening} and \emph{annotation}
#'
#'
#'
#'
#' @examples
install_db <- function()
{

  path <- paste0(system.file(package = "pato"),"/DB")
  dir.create(paste0(path))
  curl::curl_download("https://github.com/irycisBioinfo/patodb/raw/main/patodb.tar.gz", destfile = paste0(path,"/patodb.tar.gz"),quiet = F)
  untar(paste0(path,"/patodb.tar.gz"),exdir = path,compressed = "gzip")
  file.remove(paste0(path,"/patodb.tar.gz"))
}
