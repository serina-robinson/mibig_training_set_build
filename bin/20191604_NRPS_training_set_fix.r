# Install packages
pacman::p_load("DECIPHER", "plot3D", "plot3Drgl", "ade4", "muscle", "data.table",
               "phangorn", "cowplot", "RColorBrewer", "phylobase", "treeio", "ggtree", "Biostrings", "readxl", "tidyverse", "gtools", "rentrez")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Combine with the NRPS sequences
aa_dat <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
names(aa_dat) <- paste0(names(aa_dat), ".aa_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1<- word(names(aa_dat), -2, sep = "_")
subst2<- word(names(aa_dat), -3, sep = "_")
table(subst2)

# Remove all rejects
names(aa_dat)
aa_fin <- aa_dat[!grepl("reject", names(aa_dat))]
length(aa_fin) # 1072 sequences

# Fix the specificiity so that it is second
acc <- word(names(aa_fin), sep = "_", 1)
ss <- word(names(aa_fin), sep = "_", -3)
nams <- paste0(acc, "_", ss, "_", names(aa_fin))
head(nams)
names(aa_fin) <- nams
length(names(aa_fin))

# Combine the substrate specificity with the sequences from the uniprot training set
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

# 34 seqs extracted
ex34 <- readAAStringSet("output/tempfile_functional_class_34aa_extracted.faa")
names(ex34) <- sqnams_ss

# Remove those that have amino.acid as the name
ex34_noaa <- ex34[word(names(ex34), sep = "_", 2) != "amino.acid"]

# Combine the seqs
comb_sqs <- AAStringSet(c(aa_fin, ex34_noaa))
comb_dedup <- comb_sqs[!duplicated(comb_sqs)] ### How to handle duplicates - should they be removed??
length(comb_sqs)
# names(comb_dedup)
sort(table(word(names(comb_dedup), sep = "_", 2))) 

writeXStringSet(comb_dedup, "output/20191604_small_substrate_grp_for_rf_py.fasta")
