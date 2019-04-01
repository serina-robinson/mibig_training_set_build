create_colnames_15aap <-function (Alignment) 
{
  generic_list <- c("xA", "xB", "xC", "xD", "xE")
  col_names <- vector()
  Aln_length <- length((unlist(strsplit(Alignment$seq[1], split = NULL))))
  for (i in 1:Aln_length) {
    add <- sub("x", i, generic_list)
    col_names <- c(col_names, add)
  }
  return(col_names)
}