# Install packages
pacman::p_load("DECIPHER", "tidymodels", "tidyverse", "readxl")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the predictions
preds <- read_excel("data/AdenylPred-2020-04-14_35_weird_test_find_refs.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::mutate(bgcs = word(query_name, sep = "_", 1)) %>%
  dplyr::mutate(accs = word(query_name, sep = "_", -4))

preds$accs

predbgcs <- word(preds$query_name, sep = "_", 1)
predbgcs

# Read in MIBiG refs
mibig <- read_excel("data/mibig_training_set_manually_edited_20192603.xlsx") %>%
  dplyr::mutate(accs = word(acc, sep = "\\.1", 1))
# refr <- mibig[mibig$bgcs %in% predbgcs,]
mibig$acc

# Merge
comball <- preds %>%
  dplyr::left_join(., mibig, by = "accs")
dim(comball)

write_csv(comball, "output/20200414_combined_preds_data.csv")
          
# Read in Uniprotkb data
# preds$query_name
# uniprot <- read_excel("data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
# uniprot$acc

uniprot$pmid
