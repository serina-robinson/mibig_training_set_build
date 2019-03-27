## Install packages
pacman::p_load("DECIPHER", "data.table", "Biostrings", "caret", "bgafun", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
rawdat <- readAAStringSet("data/713_total_training_sqs_including_NRPS.fa")
clf <- word(names(rawdat), -1, sep = "_")
table(clf)
subst <- word(names(rawdat), -2, sep = "_")

# Ran antismash extract 34 aa extract_34_aa_serina.py script 
ex34 <- fread('../sandpuma2_serina/flat/713_aa34_signatures.fa', header = F, data.table = F)
ex34
ex34_gsub <- gsub("-", "X", ex34[662,1])
ex34_gsub
aa34 <- AAStringSet(ex34_gsub)
length(aa34)

names(aa34) <- names(rawdat)
writeXStringSet(aa34, 'data/713_aa34_gaps_removed_signatures.fa')
# readAAStringSet('data/713_aa34_signatures.fa')

# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/713_aa34_signatures.fa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)
aa <- bgafun::convert_aln_AAP(rdaln) #5 physicochemical properties

aadf <- data.frame(aa, stringsAsFactors = F)

aap <- aadf %>%
  dplyr::mutate(nms = rownames(.)) %>%
  dplyr::select(-contains("D")) %>%
  dplyr::filter(!grepl("ERROR", X1A))

rownames(aap) <- aap$nms
aap <- aap %>%
  dplyr::select(-nms)
colnames(aap) <- gsub("^X","",colnames(aap))

colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap))
numfeats <- length(colnames(aap)) # 1096
colnames(aap)


# Write to CSV file
write.csv(aap, paste0("data/713_seqs_", numfeats, "_feats_for_supervised.csv"), row.names=rownames(aap), quote = F)
# 136 sequences
