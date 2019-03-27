## Install packages
pacman::p_load('ape', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(123091)

# Read in the mibig dataset
mib <- read_excel("data/mibig_training_set_manually_edited_20192603.xlsx") %>%
  dplyr::filter(confidence > 0) %>%
  mutate(small_substrate_group = case_when(small_substrate_group == "malonate.derivatives" ~ "C2.through.C4",
                                          small_substrate_group == "salicylic.acid.derivatives" ~ "aryl.and.biaryl.derivatives",
                                           TRUE ~ small_substrate_group)) %>%
  mutate(data_source = "mibig") %>%
  dplyr::rename(protein_names = prot_names) 

# Read in the uniprot dataset
cmbnd <- read_csv("data/anl_training_set_updated_20192503_fixnams.csv") %>%
 #  dplyr::filter(!grepl("coelente", substrate)) %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  mutate(sqnams_tr = paste0(1:nrow(.), "_", org_short, "_", substrate_group, "_", functional_class)) %>%
  dplyr::mutate(likely_substrate = substrate) %>%
  dplyr::mutate(pmid = pub_med_id) %>%
  mutate(small_substrate_group = case_when(functional_class == "OTHER|HOLDOUTTEST" ~ "unknown.other",
                                           likely_substrate == "33aS4S7aS7amethyl15dioxooctahydro1Hinden4ylpropanoate" ~ "unknown.other",
                                           str_detect(likely_substrate, "carnitine|isoleucine|leuc|proli|alanine|aminolev|threonine") ~ "amino.acid",
                                           str_detect(likely_substrate, "stear|olei|olea|docosa|tetracosa|oxocholes|nonadec|arachi|eicosa|bile|verylongchain") ~ "C18.and.up.or.bile.acid",
                                           str_detect(likely_substrate, "tridec|myrist|palmit|longchain|phytodien|heptadec|long_chain") ~ "C13.through.C17",
                                           str_detect(likely_substrate, "pent|valer|hex|hept|oct|mediumchain") ~ "C5.through.C8",
                                           str_detect(likely_substrate, "non|dec|dodec|laur") ~ "C9.through.C12",
                                           str_detect(likely_substrate, "coum|cinna|vanill|ferul|succinyl|quinol") ~ "cinnamate.and.succinylbenzoate.derivatives",
                                           str_detect(likely_substrate, "luciferin") ~ "luciferin",
                                           str_detect(likely_substrate, "medianbetahydroxyacidlong") ~ "median.beta.hydroxyacid",
                                           str_detect(likely_substrate, "anthran|benzo|naphth|aryl|xanthu|quin|dihydroxybenz|salicyl|picolin|nicotin|phenylacet") ~ "aryl.and.biaryl.derivatives",
                                           str_detect(likely_substrate, "acetat|shortchain|acryl|prop|but|crot|malon|methylmal|oxala|") ~ "C2.through.C4",
                                           TRUE ~ "unknown.other")) %>%
  mutate(data_source = "uniprot_kb") %>%
  dplyr::rename(large_substrate_group = substrate_group) %>%
  bind_rows(., mib) %>%
  # mutate(small_substrate_group = case_when(small_substrate_group == "C9.through.C12" ~ "C6.through.C12"))
  dplyr::select(data_source, entry_name, acc, organism, protein_names, likely_substrate, small_substrate_group, large_substrate_group, functional_class,
  bgcs, pmid, cmpnd, pdb_id, kinetics, ec_numbers, title, sqnams_tr, aa_seq)
write_csv(data.frame(cbind(cmbnd$small_substrate_group, cmbnd$likely_substrate)), "output/uniprot_substrates_and_substrate_groups.csv")

# Fix the substrate groups a bit 
cmbnd_fix <- cmbnd %>%
  dplyr::mutate(small_substrate_group = case_when(str_detect(likely_substrate, "pentan|valer") ~ "C2.through.C5",
                                                  small_substrate_group == "C2.through.C4" ~ "C2.through.C5",
                                                  small_substrate_group == "C9.through.C12" ~ "C6.through.C12",
                                                  small_substrate_group == "C5.through.C8" ~ "C6.through.C12",
                                                  TRUE ~ small_substrate_group))
table(cmbnd_fix$small_substrate_group)
write_csv(data.frame(sort(table(cmbnd_fix$likely_substrate), decreasing = T)), "output/substrates_sorted.csv")
