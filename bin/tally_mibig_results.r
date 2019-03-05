## Install packages
pacman::p_load('ape', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set seed 
set.seed(123091)

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build/")

# Read in the MIBiG training data
rawdat <- read_excel("data/mibig_training_set_manually_edited_20190304.xlsx")
dim(rawdat)
mibig <- rawdat %>%
  dplyr::filter(confidence > 0) %>%
  mutate(sqnams = paste0(bgcs, "_", word(cmpnd, 1), "_", acc, "_", substrate_group, "_", functional_class))
mibig$sqnams <- gsub("-", "_", mibig$sqnams)
mibig$sqnams <- gsub("\\.", "_", mibig$sqnams)
dim(mibig)
table(mibig$functional_class)

mibig$likely_substrate[mibig$functional_class=="PEPTIDE"]

# Read in the UniPROT KB training data
uniprot <- read_csv("data/anl_training_set_updated_20190215_fixnams.csv") %>%
  dplyr::filter(!grepl("coelente", substrate))

uniprot$substrate_group[grep("coum", uniprot$substrate)] <- "coumarin"
uniprot$functional_class[grep("coum", uniprot$substrate)] <- "COUM"
uniprot$substrate_group[grep("naphthoate", uniprot$substrate)] <- "biaryl"
uniprot$functional_class[grep("naphthoate", uniprot$substrate)] <- "BIARYL"

# uniprot$protein_names[uniprot$functional_class == "OTHER"] # can keep class other

# Combine the MIBiG and the UNIPROT
dat <- uniprot %>%
  bind_rows(mibig, .) %>%
  dplyr::filter(!functional_class %in% c("PEPTIDE", "CAR"))
colnames(dat)

# 584 observations
table(cmbnd$functional_class)

# Pull sequences
sqs <- AAStringSet(dat$aa_seq)
names(sqs) <- dat$sqnams
names(sqs)
summary(width(sqs))

# Remove duplicates
nodups <- sqs[!duplicated(sqs)]
dups <- sqs[duplicated(sqs)]
dups
length(nodups) # 3 duplicates
nodups

# Remove sequences shorter than 400 amino acids
sqs_long <- nodups[width(nodups) >= 350] # 575 sequences
sqs_short <- nodups[width(nodups) < 350]
length(sqs_long)
# 575 sequences
writeXStringSet(sqs, "data/575_training_seqs_including_mibig.fasta")
names(sqs)

# Align sequences structurally using DECIPHER package 
aa.al <- AlignSeqs(sqs)

# Adjust alignment
aa.adj <- AdjustAlignment(aa.al)
writeXStringSet(aa.adj, "data/575_training_seqs_aligned_adjusted.fasta")

summary(width(aa.adj))
BrowseSeqs(aa.adj)
