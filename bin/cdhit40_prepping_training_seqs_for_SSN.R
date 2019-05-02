## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
rawdat <- readAAStringSet("data/669_training_set_to_extract_34aa_20192803.faa")
length(rawdat)
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid"), collapse = "|"), names(rawdat)),] # 658 observations
clf <- word(names(dat), -1, sep = "_")
subst1 <- word(names(dat), -2, sep = "_")
subst2 <- word(names(dat), -3, sep = "_")
length(dat)

# Read in the NRPS training set
aa_dat <- readAAStringSet("data/1093_sp2_full_length_names_fixed_large_grps.faa")
names(aa_dat) <- paste0(names(aa_dat), "_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1 <- word(names(aa_dat), -2, sep = "_")
subst2 <- word(names(aa_dat), -3, sep = "_")
table(subst2)
aa_rem <- aa_dat[subst2 != "reject"]

# Randomly sample 
set.seed(20190304)
ran_nums <- sample(1:length(aa_rem), 150, replace = F)
aa_dat <- aa_rem[ran_nums]

# Combine the two
comb <- AAStringSet(c(dat, aa_dat))


# Be picky about the length of the sequence
summary(width(comb))
comb2 <- comb[width(comb) < 1500]
head(comb2)
length(comb2) # 1702 
summary(width(comb2))
clf <- word(names(comb2), -1, sep = "_")
subst1<- word(names(comb2), -2, sep = "_")
subst2<- word(names(comb2), -3, sep = "_")


tabvec <- as.vector(table(clf)) 
tabvec
names(table(clf))
table(clf)

# ran_nums <- lapply(tabvec, function(x) { sample(x, size = 40, replace = F) })
# 
# ran_df <- data.frame(names(comb2), clf)
# fin_df <- ran_df %>%
#   split(.$clf)
# length(fin_df)
# fin_df[[1]]
# fin_sampled <- do.call("rbind", lapply(1:length(fin_df), function(x) {fin_df[[x]][ran_nums[[x]],]}))
# 
# finset <- comb2[as.character(fin_sampled$names.comb2.)] # 360 sequences to cluster with
# length(finset)
# # writeXStringSet(finset, "data/360_randomly_sampled_class_reps.fasta")
# summary(width(finset))
# table(word(names(finset), sep = "_", -2))

# Read in the seqs from cdhit
writeXStringSet(comb2, "data/778_functional_class_reps_from_training.fasta")
# finset <- readAAStringSet("data/360_randomly_sampled_class_reps.fasta")
finset <- readAAStringSet("data/778_functional_class_reps_from_training.fasta")
cdhit <- readAAStringSet("data/cdhit40.fa")
length(cdhit)

# Combine with the 360 sequences
comb_cdhit <- c(finset, cdhit)
length(comb_cdhit)
comb_nodups <- comb_cdhit[!duplicated(comb_cdhit)]
length(comb_nodups) # 3119
writeXStringSet(comb_nodups, "data/3319_comb_cdhit40_for_ssn.fasta")

# Write to file
# fda <- read.alignment("data/5516_comb_cdhit_for_ssn.fasta", format = "fasta")
# table(duplicated(fda$seq))

# 