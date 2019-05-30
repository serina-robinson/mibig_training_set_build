## Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190403)

# Read in the sp2 data
sp2 <- readAAStringSet('data/1093_sp2_full_length_names_fixed_large_grps.faa')
head(sp2)
table(word(names(sp2), -1, sep = "_"))

# Read in the most updated UniProt sequences and extract
uniprot <- readAAStringSet("data/669_training_set_to_extract_34aa_20192803.faa")
names(uniprot)[1:10]
table(word(names(uniprot), -3, sep = "_"))
comb2 <- comb[-grep(paste0(c("HOLDOUT", "OTHER", "amino.acid", "reject"), collapse = "|"), names(comb)),]

# Combine
comb <- c(sp2, uniprot)
length(comb)
comb2 <- comb[-grep(paste0(c("HOLDOUT", "OTHER", "amino.acid", "reject"), collapse = "|"), names(comb)),]
length(comb2)

