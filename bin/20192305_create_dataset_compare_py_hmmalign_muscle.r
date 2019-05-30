# Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "muscle", "data.table",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", 
               "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Combine with the NRPS sequences
aa_dat <- readAAStringSet("data/1093_sp2_full_length_names_fixed_large_grps.faa")
table(duplicated(aa_dat))
names(aa_dat) <- paste0(names(aa_dat), ".aa_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1<- word(names(aa_dat), -2, sep = "_")
subst2<- word(names(aa_dat), -3, sep = "_")
table(subst2)
table(clf)

# Remove all rejects
names(aa_dat)
aa_fin <- aa_dat[!grepl("reject", names(aa_dat))]
length(aa_fin) # 1072 sequences

# Fix the specificiity so that it is second
acc <- word(names(aa_fin), sep = "_", 1)
ss <- word(names(aa_fin), sep = "_", -3)
nams <- paste0(acc, "_", ss, "_NRPS_", names(aa_fin))
names(aa_fin) <- nams
table(word(names(aa_fin), sep = "_", 1))
table(word(names(aa_fin), sep = "_", 2))
table(word(names(aa_fin), sep = "_", 3))


# Combine the substrate specificity with the sequences from the uniprot training set
rawdat <- read_excel("../mibig_training_set_build_test/data/combined_adenylate_forming_training_set_for_db_20191404.xlsx") %>%
  dplyr::filter(!functional_class %in% c("OTHER", "CAR", "NRPS")) %>%
  dplyr::filter(small_substrate_group != "unknown.other") %>%
  dplyr::mutate(org_clned = paste0(word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "))) %>%
  dplyr::mutate(acc_clned = gsub("_", "", acc))
table(rawdat$functional_class, rawdat$small_substrate_group)
rawdat$org_clned
rawdat$acc_clned

sqs <- AAStringSet(rawdat$aa_seq)
length(sqs)
# sqnams_ss <- gsub(" ", "_", paste0(rawdat$acc, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group))
# sqnams_fc <- gsub(" ", "_", paste0(rawdat$acc, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group, "_", rawdat$functional_class))
sqnams_fc <- gsub(" ", "_", paste0(rawdat$acc_clned, "_", rawdat$functional_class, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$small_substrate_group))
sqnams_ss <- gsub(" ", "_", paste0(rawdat$acc_clned, "_",  rawdat$small_substrate_group, "_", rawdat$functional_class, "_", rawdat$org_clned, "_", rawdat$small_substrate_group, "_", rawdat$functional_class))

names(sqs) <- sqnams_ss
table(word(names(sqs), 3, sep = "_"))
head(names(sqs))
length(sqs)

# Fix names
writeXStringSet(sqs, "~/Documents/Wageningen_UR/github/amplicon_pred/data/625_uniprot_seqs_subspec_nams_no_NRPS.faa")


# Combine the seqs
comb_sqs <- AAStringSet(c(aa_fin, sqs))
length(comb_sqs)
summary(width(comb_sqs))
comb_dedup <- comb_sqs[!duplicated(comb_sqs)] ### How to handle duplicates - should they be removed??
length(comb_dedup)
sort(table(word(names(comb_sqs), sep = "_", 2))) 

writeXStringSet(comb_dedup, "output/20192305_1693_small_substrate_grp_for_rf_duplicates_removed.fasta")
writeXStringSet(comb_dedup, "~/Documents/Wageningen_UR/github/amplicon_pred/data/20192305_1693_small_substrate_grp_for_rf_duplicates_removed.fasta")
