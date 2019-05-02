## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190304)

# Read in the training set
# origdat <- readAAStringSet("data/676_training_set_to_extract_34aa_20192703.faa")
rawdat <- readAAStringSet("data/669_training_set_34aa_extracted_20192703.faa")
length(rawdat)
clf <- word(names(rawdat), -1, sep = "_")
subst1<- word(names(rawdat), -2, sep = "_")
subst2<- word(names(rawdat), -3, sep = "_")

# Read in the NRPS training set
aa_dat <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
ran_nums <- sample(1:length(aa_dat), 150, replace = F)
aa_dat <- aa_dat[ran_nums]

names(aa_dat) <- paste0(names(aa_dat), "_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1<- word(names(aa_dat), -2, sep = "_")
subst2<- word(names(aa_dat), -3, sep = "_")
table(subst2)
aa_rem <- aa_dat[subst2 != "reject"]

# Combine the two
comb <- AAStringSet(c(rawdat, aa_rem))
length(comb)
# writeXStringSet(comb, "data/817_aa34_signatures_20190104.fa")

# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/817_aa34_signatures_20190104.fa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)

source("src/convert_aln_15aap.r")
aa <- convert_aln_15aap(rdaln) #15 physicochemical properties
aadf <- data.frame(aa, stringsAsFactors = F)

aap <- aadf %>%
  dplyr::mutate(nms = rownames(.)) %>%
  # dplyr::select(-contains("D")) %>%
  # dplyr::filter(!grepl("ERROR", X1_WOLS870101)) %>%
  dplyr::filter(word(nms, -3, sep = "_") != "amino.acid")

rownames(aap) <- aap$nms

# table(word(aap$nms, -3, sep = "_"))

aap <- aap %>%
  dplyr::select(-nms)
colnames(aap) <- gsub("^X","",colnames(aap))

# colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap))
numfeats <- length(colnames(aap)) # 1096
colnames(aap)
dim(aap)

# Write to CSV file
write.csv(aap, paste0("data/797_seqs_510_feats_for_supervised_centered_scaled_20191104.csv"), row.names=rownames(aap), quote = F)

# Read in the data
# rawdat <- read_csv("data/797_seqs_510_feats_for_supervised_centered_scaled_20191104.csv")
