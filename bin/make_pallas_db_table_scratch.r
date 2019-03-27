## Install packages
pacman::p_load('ape', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set seed 
set.seed(123091)

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the MIBiG training data
rawdat <- read_excel("data/mibig_training_set_manually_edited_20190304.xlsx")

mibig <- rawdat %>%
  dplyr::filter(confidence > 0) %>%
  mutate(substrate_group = str_replace_all(substrate_group, "coumarin", "aryl")) %>%
  mutate(substrate_group = str_replace_all(substrate_group, "biaryl", "aryl")) %>%
  mutate(substrate_group_tr = str_replace_all(substrate_group, "_", "")) %>%
  mutate(functional_class = str_replace_all(functional_class, "BIARYL", "ARYL")) %>%
  mutate(functional_class = str_replace_all(functional_class, "COUM", "ARYL")) %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  mutate(sqnams_tr = paste0(bgcs, "_", word(cmpnd, 1, sep = "_"), "_", acc, "_", substrate_group_tr, "_", functional_class)) %>%
  mutate(sqnams_tr = str_replace_all(sqnams_tr, "-", "_")) %>%
  mutate(sqnams_tr = str_replace_all(sqnams_tr, "\\.", "_"))

# Read in the UniPROT KB training data
uniprot <- read_csv("data/anl_training_set_updated_20190215_fixnams.csv") %>%
  dplyr::filter(!grepl("coelente", substrate)) %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  mutate(sqnams_tr = paste0(1:nrow(.), "_", org_short, "_", substrate_group, "_", functional_class))
uniprot$sqnams_tr

# Combine the MIBiG and the UNIPROT
# uniprot <- read_excel('data/anl_training_set_updated_20190215.xlsx')
dat <- uniprot %>%
  dplyr::mutate(likely_substrate = substrate) %>%
  dplyr::mutate(pmid = pub_med_id) %>%
  bind_rows(., mibig) %>%
  dplyr::select(entry_name, organism, protein_names, likely_substrate, substrate_group, functional_class,
                bgcs, pmid, cmpnd, pdb_id, kinetics, ec_numbers, title, sqnams_tr, aa_seq) #%>%

dat_seqs <- AAStringSet(dat$aa_seq)
names(dat_seqs) <- gsub(" ", "_", dat$sqnams_tr)
table(word(names(dat_seqs), sep = "_", -1))

writeXStringSet(dat_seqs, '~/Documents/Wageningen_UR/github/antismash_clustering/data/613_training_set_seqs.fa')

# write_csv(dat, "data/20190318_full_anl_training_set.csv")

# write_csv(dat, "../pallas/shiny_app/data/20190318_full_anl_training_set.csv")

#  dplyr::filter(!functional_class %in% c("PEPTIDE", "CAR"))
colnames(dat)
table(dat$substrate)
# 584 observations
table(dat$functional_class)

# Pull sequences
sqs <- AAStringSet(dat$aa_seq)
names(sqs) <- dat$sqnams_tr
names(sqs)
summary(width(sqs))

# Remove duplicates
nodups <- sqs[!duplicated(sqs)]
dups <- sqs[duplicated(sqs)]
names(dups)
length(nodups) # 3 duplicates
nodups

# Remove sequences shorter than 400 amino acids
sqs_long <- nodups[width(nodups) >= 350] # 575 sequences
sqs_short <- nodups[width(nodups) < 350]
length(sqs_long) # 602

# 575 sequences
writeXStringSet(sqs, "data/602_training_seqs_including_CAR_PEPTIDE.fasta")
names(sqs)

