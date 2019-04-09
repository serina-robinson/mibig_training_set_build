## Install packages
pacman::p_load("caret", "readxl", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

colnames(apred)
# Read in the adenylpred
apred <- read_csv("data/AdenylPred-2019-04-09.csv") %>%
  janitor::clean_names() %>%
  mutate(truth_subspec = word(query_name, sep = "_", -3)) %>%
  mutate(truth_class = word(query_name, sep = "_", -1))

apred
write_csv(apred, "data/adenylpred_output_20190409.csv")
