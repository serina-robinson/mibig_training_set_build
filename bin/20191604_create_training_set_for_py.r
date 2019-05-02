# Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "muscle", "data.table",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the full length sequences
rawdat <- read_excel("../mibig_training_set_build_test/data/combined_adenylate_forming_training_set_for_db_20191404.xlsx") %>%
  dplyr::filter(!functional_class %in% c("OTHER", "CAR")) %>%
  dplyr::filter(small_substrate_group != "unknown.other") %>%
  dplyr::mutate(org_clned = paste0(word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "))) %>%
  dplyr::mutate(acc_clned = gsub("_", "", acc))
table(rawdat$functional_class, rawdat$small_substrate_group)
rawdat$org_clned
rawdat$acc_clned

sqs <- AAStringSet(rawdat$aa_seq)
# sqnams_ss <- gsub(" ", "_", paste0(rawdat$acc, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group))
# sqnams_fc <- gsub(" ", "_", paste0(rawdat$acc, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group, "_", rawdat$functional_class))
sqnams_fc <- gsub(" ", "_", paste0(rawdat$acc_clned, "_", rawdat$functional_class, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group))
sqnams_ss <- gsub(" ", "_", paste0(rawdat$acc_clned, "_",  rawdat$small_substrate_group, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$functional_class))


names(sqs) <- sqnams_fc
head(names(sqs))
# subst2<- word(names(sqs), -3, sep = "_")


# Extract the 34 amino acids
source("src/extract_34_aa.r")
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa(query_fils[x]) })


# Combine with the NRPS sequences
aa_dat <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
names(aa_dat) <- paste0(names(aa_dat), ".aa_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1<- word(names(aa_dat), -2, sep = "_")
subst2<- word(names(aa_dat), -3, sep = "_")
table(subst2)
aa_rem <- aa_dat[subst2 != "reject.aa"]
names(aa_rem)

# 