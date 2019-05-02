sqs <- readAAStringSet(input$file1$datapath)
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa(query_fils[x]) })
extract_34_df <- data.frame(matrix(unlist(extract_34_list), nrow = length(sqs), byrow=T), 
                            stringsAsFactors=FALSE)
colnames(extract_34_df) <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
