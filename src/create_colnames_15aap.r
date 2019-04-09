create_colnames_15aap <-function (Alignment) 
{
  # generic_list <- c("xA", "xB", "xC", "xD", "xE")
  tmp <- read_csv("data/15_aa_physicochemical_properties.csv")
  col_names <- paste0(colnames(tmp)[3:17], "_", rep(1:34, each = 15))
  # col_names <- vector()
  # Aln_length <- length((unlist(strsplit(Alignment$seq[1], split = NULL))))
  # for (i in 1:Aln_length) {
  #   add <- sub("x", i, generic_list)
  #   col_names <- c(col_names, add)
  # }
  return(col_names)
}