## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the test set IDs
ids <- fread("output/test_set_seqs_for_tool_prediction.csv", sep = ",", data.table = F)[,1]
ids

# Search through files to extract test set sequences
fa1 <- readAAStringSet("data/669_training_set_to_extract_34aa_20192803.faa")
fa1_test <- fa1[names(fa1) %in% ids]
length(fa1_test)
writeXStringSet(fa1_test, "output/141_training_set_seqs.faa")
# 