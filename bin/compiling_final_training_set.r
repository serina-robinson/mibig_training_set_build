## Install packages
pacman::p_load('ape', 'tidyverse', 'data.table', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in mibig organisms
mib_orgs <- read_excel('data/20190327_mibig_organisms_bgcs.xlsx') %>%
  janitor::clean_names() %>%
  mutate(org_short = paste0(word(organism, sep = " ", 1), "_", word(organism, sep = " ", 2))) %>%
  dplyr::select(-organism) %>%
  dplyr::rename(bgcs = mi_bi_g_accession,
                organism = org_short) %>%
  dplyr::select(bgcs, organism)

# Read in the mibig dataset
mib <- read_excel("data/mibig_training_set_manually_edited_20192603.xlsx") %>%
  dplyr::filter(confidence > 0) %>%
  left_join(., mib_orgs) %>%
  mutate(small_substrate_group = case_when(small_substrate_group == "malonate.derivatives" ~ "C2.through.C4",
                                          small_substrate_group == "salicylic.acid.derivatives" ~ "aryl.and.biaryl.derivatives",
                                           TRUE ~ small_substrate_group)) %>%
  mutate(data_source = "mibig") %>%
  dplyr::rename(protein_names = prot_names) %>%
  dplyr::mutate(old_sq_nams = sqnams_tr) %>%
  dplyr::mutate(small_substrate_group = case_when(str_detect(likely_substrate, "pentan|valer") ~ "C2.through.C5",
                                                  small_substrate_group == "C2.through.C4" ~ "C2.through.C5",
                                                  small_substrate_group == "C9.through.C12" ~ "C6.through.C12",
                                                  small_substrate_group == "C5.through.C8" ~ "C6.through.C12",
                                                  TRUE ~ small_substrate_group)) %>%
  dplyr::mutate(large_substrate_group = case_when(small_substrate_group == "C2.through.C5" ~ "shortchain",
                                                  small_substrate_group == "C6.through.C12" ~ "mediumchain",
                                                  small_substrate_group == "C13.through.C17" ~ "longchain",
                                                  small_substrate_group == "C18.and.up.or.bild.acid" ~ "verylongchainbile",
                                                  TRUE ~ large_substrate_group)) %>%
  dplyr::mutate(sqnams_tr = paste0(bgcs, "_", acc, "_", word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "), "_",
                                   small_substrate_group, "_", large_substrate_group, "_", functional_class))

# Read in the uniprot dataset
cmbnd <- read_excel("data/anl_training_set_updated_20192703.xlsx") %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  #mutate(sqnams_tr = paste0(1:nrow(.), "_", org_short, "_", substrate_group, "_", functional_class)) %>%
  dplyr::mutate(likely_substrate = substrate_condensed) %>%
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
  dplyr::mutate(old_sq_nams = sqnams_tr) %>%
  dplyr::mutate(entry = case_when(functional_class == "BLS" ~ word(protein_names, 1, sep = " "),
                                  TRUE ~ entry)) %>%
  dplyr::mutate(gene_names = case_when(functional_class == "BLS" ~ "polyolefin biosynthesis oleC",
                                                TRUE ~ gene_names)) %>%
  dplyr::mutate(entry_name = case_when(functional_class == "BLS" ~ word(protein_names, 1, sep = " "),
                                       TRUE ~ entry_name)) %>%
  dplyr::rename(acc = entry_name)


inds <- which(cmbnd$acc == "NA")
cmbnd$acc[inds] <- word(cmbnd$protein_names[inds], sep = " ", 1)
cmbnd$functional_class[cmbnd$functional_class == "PEPTIDE"] <- "NRPS"
cmbnd$functional_class[cmbnd$functional_class == "VLACS_BILE"] <- "VLACSBILE"
cmbnd$functional_class[cmbnd$functional_class == "FAT"] <- "VLACSBILE"
cmbnd$small_substrate_group[cmbnd$functional_class == "VLACSBILE"] <- "C18.and.up.or.bile.acid"
cmbnd$large_substrate_group[cmbnd$functional_class == "VLACSBILE"] <- "verylongchainbile"

vlacs <- c("92_Drosophila_melanogaster_longchain_LACS", "153_Mycobacterium_tuberculosis_longchain_FAAL", "138_Sus_scrofa_longchain_LACS",
           "128_Mycobacterium_tuberculosis_longchain_VLACSBILE", "155_Mycobacterium_tuberculosis_longchain_LACS")



cmbnd$likely_substrate[cmbnd$old_sq_nams %in% vlacs] <- "C18.and.up.or.bile.acid"
cmbnd$substrate[cmbnd$old_sq_nams %in% vlacs] <- "C18.and.up.or.bile.acid"
cmbnd$small_substrate_group[cmbnd$old_sq_nams  %in% vlacs] <- "C18.and.up.or.bile.acid"
cmbnd$large_substrate_group[cmbnd$old_sq_nams  %in% vlacs] <- "verylongchainbile"
cmbnd$functional_class[cmbnd$old_sq_nams  %in% vlacs] <- "VLACSBILE"


# cmbnd$likely_substrate[cmbnd$old_sq_nams == "138_Sus_scrofa_longchain_LACS"] <- "C18.and.up.or.bile.acid"
# cmbnd$small_substrate_group[cmbnd$old_sq_nams == "138_Sus_scrofa_longchain_LACS"] <- "C18.and.up.or.bile.acid"
# cmbnd$functional_class[cmbnd$old_sq_nams == "138_Sus_scrofa_longchain_LACS"] <- "VLACSBILE"


# write_csv(data.frame(cbind(cmbnd$small_substrate_group, cmbnd$likely_substrate)), "output/uniprot_substrates_and_substrate_groups.csv")

# Fix the substrate groups a bit 
cmbnd_fix <- cmbnd %>%
  dplyr::mutate(small_substrate_group = case_when(str_detect(likely_substrate, "pentan|valer") ~ "C2.through.C5",
                                                  small_substrate_group == "C2.through.C4" ~ "C2.through.C5",
                                                  small_substrate_group == "C9.through.C12" ~ "C6.through.C12",
                                                  small_substrate_group == "C5.through.C8" ~ "C6.through.C12",
                                                  TRUE ~ small_substrate_group)) %>%
  dplyr::mutate(large_substrate_group = case_when(small_substrate_group == "C2.through.C5" ~ "shortchain",
                                                  small_substrate_group == "C6.through.C12" ~ "mediumchain",
                                                  small_substrate_group == "C13.through.C17" ~ "longchain",
                                                  small_substrate_group == "C18.and.up.or.bile.acid" ~ "verylongchainbile",
                                                  TRUE ~ large_substrate_group)) %>%
  dplyr::mutate(sqnams_tr = paste0(acc, "_", word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "), "_",
                                   small_substrate_group, "_", large_substrate_group, "_", functional_class)) %>%
  bind_rows(., mib) %>%
  dplyr::filter(functional_class != "HOLDOUTTEST") %>%
  # mutate(small_substrate_group = case_when(small_substrate_group == "C9.through.C12" ~ "C6.through.C12"))
  dplyr::select(data_source, acc, acc, organism, protein_names, likely_substrate, substrate, small_substrate_group, large_substrate_group, functional_class,
                bgcs, pmid, cmpnd, pdb_id, kinetics, ec_numbers, title, sqnams_tr, aa_seq)


# write_csv(data.frame(sort(table(cmbnd_fix$likely_substrate), decreasing = T)), "output/substrates_sorted.csv")

# Need to add more coumarate-CoA ligases to the dataset
coums <- fread("data/uniprot-coumar_+CoA+ligase+existence__Evidence+at+transcript+level+[2]%2--.tab", sep = "\t", data.table = F)
coums_nodup <- coums[!duplicated(coums$Sequence),] %>%
  dplyr::filter(Length > 350)

coums_aa <- AAStringSet(coums_nodup$Sequence)
names(coums_aa) <- gsub(" ", "_", paste0(coums_nodup$`Entry name`, "_", coums_nodup$Organism))
# writeXStringSet(coums_aa, "data/251_coums_for_clustering.fa")

coums_cdhit <- readAAStringSet('data/coums_condensed_for_clustering_0.7.fa')
tofind <- paste0(word(names(coums_cdhit), sep = "_", 1), "_", word(names(coums_cdhit), sep = "_", 2)) 
seqs_include <- coums_nodup[coums_nodup$`Entry name` %in% tofind,]
seqs_add <- seqs_include %>%
  janitor::clean_names() %>%
  mutate(data_source = "uniprot_kb",
         likely_substrate = "4coumarate",
         substrate = "4-coumarate",
         small_substrate_group = "cinnamate.and.succinylbenzoate.derivatives",
         large_substrate_group = "aryl",
         functional_class = "ARYL",
         sqnams_tr = paste0(entry_name, "_", word(organism, 1, sep = " "), "_", word(organism, 2, sep = " "), "_", small_substrate_group, "_", large_substrate_group, "_", functional_class)
         ) %>%
  dplyr::rename(pmid = pub_med_id,
                aa_seq = sequence) %>%
  dplyr::select(data_source, entry_name, organism, protein_names, likely_substrate, substrate, small_substrate_group, large_substrate_group, functional_class,
                pmid, kinetics, sqnams_tr, aa_seq) %>%
  dplyr::rename(acc = entry_name)


cmbnd_coums <- cmbnd_fix %>%
  bind_rows(seqs_add)

table(cmbnd_coums$small_substrate_group)
table(cmbnd_coums$large_substrate_group)
table(cmbnd_coums$functional_class)
# cmbnd_coums[cmbnd_coums$functional_class == "VLACS_BILE",]

cmbnd_coums$organism <- gsub("_", " ", cmbnd_coums$organism)
cmbnd_coums$sqnams_tr <- gsub(" ", "_", cmbnd_coums$sqnams_tr)
dim(cmbnd_coums)
# write_csv(cmbnd_coums, "data/combined_adenylate_forming_training_set_20192703.csv")

