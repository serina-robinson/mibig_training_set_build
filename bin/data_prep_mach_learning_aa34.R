## Install packages
pacman::p_load("DECIPHER", "data.table", "Biostrings", "caret", "bgafun", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
# origdat <- readAAStringSet("data/676_training_set_to_extract_34aa_20192703.faa")
rawdat <- readAAStringSet("data/669_training_set_34aa_extracted_20192703.faa")
length(rawdat)
clf <- word(names(rawdat), -1, sep = "_")
subst1<- word(names(rawdat), -2, sep = "_")
subst2<- word(names(rawdat), -3, sep = "_")

# Read in the NRPS training set
aa_dat <- readAAStringSet("data/adom_signatures_one_per_class_34aa_20190328.faa")
names(aa_dat) <- paste0(names(aa_dat), "_amino.acid_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1<- word(names(aa_dat), -2, sep = "_")
subst2<- word(names(aa_dat), -3, sep = "_")

# Combine the two
comb <- AAStringSet(c(rawdat, aa_dat))
length(comb)
writeXStringSet(comb, "data/793_aa34_signatures_20190329.fa")

# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/793_aa34_signatures_20190329.fa', format = "fasta")
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
dim(aap)

dim(aap)
# Write to CSV file
write.csv(aap, paste0("data/793_seqs_136_feats_for_supervised_20190329.csv"), row.names=rownames(aap), quote = F)
