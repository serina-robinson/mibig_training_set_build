## Install packages
pacman::p_load('ape', 'tidyverse', 'data.table', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
adoms <- readAAStringSet("../sandpuma2_bitbucket/flat/sp2.adomains.faa")

# A-domain classes 
subst <- word(names(adoms), sep = "\t", 2) 
acc <- word(names(adoms), sep = "\t", 3) 

