# Install packages
pacman::p_load("DECIPHER", "tidymodels", "tidyverse", "readxl", "data.table")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the HMM output table
hmm <- fread("output/AMP_binding_hits_mibig_since_2019.txt", skip = 3, fill = T, data.table = F)
hmm2 <- hmm %>%
  group_by(V1) %>%
  slice(1) %>%
  dplyr::filter(!grepl("ribosomal|polyketide|NRPS|#", V1)) %>%
  dplyr::mutate(bgc = word(V1, sep = "\\|", 1)) %>%
  dplyr::mutate(accs = word(V1, sep = "\\|", 5))
  

# Read in metadata
metadat <- read_excel("output/new_BGCs_since_2019_table.xlsx", col_names = F)
colnames(metadat) <- c("bgc", "product", "organism")

# Combine
comb <- hmm2 %>%
  dplyr::left_join(., metadat, by = "bgc") %>%
  dplyr::select(V1, bgc, accs, product, organism)
write_csv(comb, "output/new_BGCs_non_NRPS_with_organisms.csv")

# Find these sequences and predict with AdenylPred
tofind <- comb$accs
tofind
fast <- readAAStringSet("data/MIBiG_accessions_since_2019.fasta")
finr <- fast[grepl(paste0(tofind, collapse = "|"), names(fast))]
width(finr)
writeXStringSet(finr, "output/full_length_MIBiG_new_BGCs_Adomains_not_extracted.fasta")

# Compare with training set
trset <- read_csv("data/1553_training_sqs_with_loop_extracted.csv") %>%
  select(X1) %>%
  dplyr::mutate(bgcraw = word(word(X1, sep = "BGC", 2), sep = "_", 1)) %>%
  dplyr::mutate(bgc = paste0("BGC", bgcraw))
intersect(trset$bgc, comb$bgc) # luckily no overlap with training set

# Read in AdenylPred predictions and compare head to head
rpred <- read_csv("data/AdenylPred-2020-04-14_predictions_for_new_MIBiG_Adomains") %>%
  janitor::clean_names() %>%
  dplyr::mutate(bgc = word(query_name, sep = "\\|", 1)) %>%
  dplyr::mutate(accs = word(query_name, sep = "\\|", 5))

pypred <- fread("output/adenylpred_python_commandline_predicitons_new_BGCs_Adomains_extracted.txt", data.table = F) %>%
  janitor::clean_names() 
colnames(pypred) <- paste0("py_", colnames(pypred))

pydat <- pypred %>%
  dplyr::mutate(bgc = word(py_query_name, sep = "\\|", 1)) %>%
  dplyr::mutate(accs = word(py_query_name, sep = "\\|", 5))


comball <- comb %>%
  ungroup() %>%
  dplyr::left_join(., rpred, by = "accs") %>%
  dplyr::left_join(., pydat, by = "accs") %>%
  dplyr::select(-bgc, -bgc.y, -V1) %>%
  dplyr::rename(bgc = bgc.x)
colnames(comball)
dim(comball) # 146 by 15

write_csv(comball, "output/20200414_predictions_for_new_mibig_full.csv")
