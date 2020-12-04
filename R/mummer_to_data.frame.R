mummer_to_data.frame <- function(file)
{

  coords = fread(file) %>% as_tibble()

  colnames(coords) = c("StartRef","EndRef","StartQuery",
                       "EndQuery","LenRef","LenEnd",
                       "Identity",
                       "RefLength","QueryLength",
                       "CovRef","CovQuery","Ref","Query")
  tmp = data.frame()
  for(i in 1:nrow(coords))
  {
    tmp = bind_rows(tmp, data.frame(pos = seq(coords$StartRef[i]:coords$EndRef[i]), contig = coords$Ref[i]))
  }
  return(tmp)
}

