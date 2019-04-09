convert_aln_15aap <- function(Alignment){
  ## Convert Alignment into matrix using 15 AAindex amino acid properties
  source("src/convert_seq_15aap.r")
  source("src/create_colnames_15aap.r")
  seq_count<-length(Alignment$seq)
  z<-list()
  
  for(i in 1:seq_count){
    x<- (unlist(strsplit(unlist(Alignment$seq[i]),split=NULL)))
    y<-convert_seq_15aap(x)
    z[[i]]<-y
  }
  
  z<-(as.data.frame(z))
  z<-t(z)
  rownames(z)<-Alignment$nam
  colnames(z)<-create_colnames_15aap(Alignment)
  
  return(z)
}