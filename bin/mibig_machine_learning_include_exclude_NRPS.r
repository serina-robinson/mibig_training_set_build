## Install packages
pacman::p_load('ape', 'tidymodels', 'data.table', 'caret', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')
library(devtools)
source_gist("https://gist.github.com/mrdwab/6424112")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build/")

# Read in the 575 training set seqs
training_sqs <- readAAStringSet('data/575_training_seqs_including_mibig.fasta')
names(training_sqs)

# wget https://bitbucket.org/chevrm/sandpuma2/raw/d162e561753a43834c242b9ed6dcc66849a50046/flat/sp2.adomainsD.faa
# data/sp2.adomainsD.faa

# Pull in Marc's training data
nrps <- readAAStringSet("data/sp2.adomainsD.faa")
monomers <- word(names(nrps), 2, sep = "\t")
names(nrps) <- paste0(monomers, "_", names(nrps), "_NRPS")
names(nrps) <- gsub("\\\t", "_", names(nrps))


mondf <- data.frame(cbind(monomers, names(nrps)), stringsAsFactors = F) %>%
  dplyr::group_by(monomers) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
mondf$monomers

write.csv(matrix(mondf$monomers, ncol = 13), 'data/monomers_table.csv')
table(mondf$monomers) # take one from each monomer
nrps_rep <- nrps[names(nrps) %in% mondf$V2]
length(nrps_rep)
width(nrps_rep)

# Search for duplicates with training sqs
cmbnd <- AAStringSet(c(training_sqs, nrps_rep))
nodups <- cmbnd[!duplicated(cmbnd)]

# Align sequences structurally using DECIPHER package 
aa.al <- AlignSeqs(nodups)

# Adjust alignment
aa.adj <- AdjustAlignment(aa.al)
length(aa.adj)
writeXStringSet(aa.adj, "data/710_training_seqs_aligned_adjusted.fasta")

# Trim all columns that have a gap using the trimAl tool
# http://trimal.cgenomics.org/
# ./trimal -in ~/Documents/Wageningen_UR/github/mibig_training_set_build/data/710_training_seqs_aligned_adjusted.fasta -gt 0.9 -out ~/Documents/Wageningen_UR/github/mibig_training_set_build/data/710_trainign_seqs_trimal_90.fasta

# Read in the trimmed alignment
trimaln <- readAAStringSet("data/710_training_seqs_trimal_90.fasta")
width(trimaln) # 398
# BrowseSeqs(trimaln)

rdaln <- read.alignment(file = "data/710_training_seqs_trimal_90.fasta", format = "fasta")
rdaln$seq <- toupper(rdaln$seq)
aa <- bgafun::convert_aln_AAP(rdaln) #5 physicochemical properties

# Convert to AAP
aap <- data.frame(aa) %>%
  select(-contains("D"))
colnames(aap) <- gsub("^X","",colnames(aap))
colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap))

# A: Polarity index
# B: Secondary structure
# C: Molecular size
# D: Number of codons # Remove D's!!!
# E: Electrostatic charge

# Remove 100% sequences or those that do not vary much
# sds <- apply(aap[,2:ncol(aap)], 2, sd)
# rm <- aap[,sds > 0.25]
# dim(rm)

# Write to csv file
write.csv(aap, "data/710_seqs_1592_feats_for_supervised.csv", row.names=rownames(aap), quote = F)



