# Install packages
# pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "data.table",
#                "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")
# 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
pacman::p_load("tidyverse", "Biostrings")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190304)

# Read in the dataset
rawdat <- read_csv("output/655_uniprot_train_loop_extracted.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")

# Combine with the NRPS sequences
aa_dat <- readAAStringSet("data/sp2.adomains.faa")
head(names(aa_dat))
aa_grps <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
names(aa_dat) <- paste0("NRPS_NRPS_", names(aa_grps))
sqs <- aa_dat

# Extract for all sequences
source("src/extract_34_aa_loop.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa_loop(query_fils[x]) })
extract_34_df_NRPS <- data.frame(matrix(unlist(extract_34_list), nrow = length(extract_34_list), byrow=T), stringsAsFactors=FALSE)
colnames(extract_34_df_NRPS)[1:510] <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
rownames(extract_34_df_NRPS) <- names(sqs)
write.csv(extract_34_df_NRPS, "output/1093_NRPS_train_loop_extracted.csv", quote = F, row.names = T)
dim(extract_34_df_NRPS) # 1093 x 585

# Combine the NRPS and the Uniprot
comb_df <- read_csv("output/655_uniprot_train_loop_extracted.csv") %>%
  column_to_rownames("X1")

comb_df2 <- read_csv("output/1093_NRPS_train_loop_extracted.csv") %>%
  column_to_rownames("X1")
total_dat <- data.frame(rbind(comb_df, comb_df2), stringsAsFactors = F)
rownames(total_dat)

table(duplicated(total_dat)) # 102 duplicates
total_dedup <- total_dat[!duplicated(total_dat),]
dim(total_dedup)
table(word(rownames(total_dedup), -1, sep = "_"))
dim(total_dedup)

write.csv(total_dedup, "data/1553_training_sqs_with_loop_extracted.csv", quote = F, row.names = T) 


