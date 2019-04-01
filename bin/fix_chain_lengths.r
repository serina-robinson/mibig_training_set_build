## Install packages
pacman::p_load("caret", "readxl", "Biostrings", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the most recent training set
trdat <- read_excel("data/combined_adenylate_forming_training_set_20192803.xlsx")
table(trdat$substrate == "NA")
cmbdat <- trdat %>%
  dplyr::mutate(small_substrate_group = case_when(grepl("laur|dodec|tetradeca", likely_substrate) ~ "C12.through.C17",
                                                  small_substrate_group == "C13.through.C17" ~ "C12.through.C17",
                                                  small_substrate_group == "C6.through.C12" ~ "C6.through.C11",
                                                  TRUE ~ small_substrate_group)) %>%
  dplyr::mutate(large_substrate_group = case_when(small_substrate_group == "C2.through.C5" ~ "shortchain",
                                                  small_substrate_group == "C6.through.C11" ~ "mediumchain",
                                                  small_substrate_group == "C12.through.C17" ~ "longchain",
                                                  small_substrate_group == "C18.and.up.or.bile.acid" ~ "verylongchainbile",
                                                  TRUE ~ large_substrate_group)) %>%
  dplyr::mutate(functional_class = case_when(functional_class == "FAAL" ~ "FAAL",
                                             small_substrate_group == "C6.through.C11" ~ "MACS",
                                             small_substrate_group == "C12.through.C17" ~ "LACS",
                                             TRUE ~ functional_class)) %>%
  dplyr::mutate(sqnams_tr = paste0(acc, "_", word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "), "_",
                                   small_substrate_group, "_", large_substrate_group, "_", functional_class))


cmbdat <- cmbdat %>%
  dplyr::filter(!grepl("E7EPS5_HUMAN|B7Z225_HUMAN|E7FCY5_DANRE|ACSL1_HUMAN", acc))
cmbdat$likely_substrate[cmbdat$likely_substrate == "decanote"] <- "decanoate"
cmbdat$substrate[cmbdat$likely_substrate == "decanote"] <- "decanoate"

# write_csv(cmbdat[grep("deca", cmbdat$likely_substrate),], "data/check_decanoates.csv")
write_csv(cmbdat, "data/combined_adenylate_forming_training_set_20192903.csv")

trn <- read_csv('data/combined_adenylate_forming_training_set_20192903.csv')

table(word(trn$sqnams_tr, -1, sep = "_"))
table(word(trn$sqnams_tr, -2, sep = "_"))
table(word(trn$sqnams_tr, -3, sep = "_"))

# Make AAStringSet
trn_aa <- AAStringSet(trn$aa_seq)
names(trn_aa) <- trn$sqnams_tr
trn_nodups <- trn_aa[!duplicated(trn_aa)]
length(trn_nodups)


writeXStringSet(trn_aa, "~/Documents/Wageningen_UR/github/sandpuma2_serina/flat/669_training_set_to_extract_34aa_20192803.faa")

### Noise to remove
# Danio rerio sequences?
# E7EPS5_HUMAN
# 
