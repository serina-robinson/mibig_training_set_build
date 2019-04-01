## Install packages
pacman::p_load('ape', 'tidyverse', 'data.table', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
trn <- read_csv('data/combined_adenylate_forming_training_set_20192703.csv')

table(word(trn$sqnams_tr, -1, sep = "_"))
table(word(trn$sqnams_tr, -2, sep = "_"))
table(word(trn$sqnams_tr, -3, sep = "_"))

# Make AAStringSet
trn_aa <- AAStringSet(trn$aa_seq)
names(trn_aa) <- trn$sqnams_tr
writeXStringSet(trn_aa, "data/676_training_set_to_extract_34aa_20192703.faa")
