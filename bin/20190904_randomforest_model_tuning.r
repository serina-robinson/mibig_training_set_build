## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/1742_aa34_signatures_20190104.fa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)

source("src/convert_aln_15aap.r")
aa <- convert_aln_15aap(rdaln) #5 physicochemical properties
aadf <- data.frame(aa,stringsAsFactors = F)