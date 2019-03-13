## Install packages
pacman::p_load('ape', 'tidymodels', 'stringr', 'data.table', 'caret', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')
# library(devtools)
# source_gist("https://gist.github.com/mrdwab/6424112")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the 575 training set seqs
training_sqs <- readAAStringSet('data/602_training_seqs_including_CAR_PEPTIDE.fasta')
table(word(names(training_sqs), -2, sep = "_"))

# wget https://bitbucket.org/chevrm/sandpuma2/raw/d162e561753a43834c242b9ed6dcc66849a50046/flat/sp2.adomainsD.faa
# data/sp2.adomainsD.faa

# Read in the JGI sequences and append
jgi <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl("OleC", role))
jgi_sqs <- AAStringSet(jgi$seq)
names(jgi_sqs) <- paste0(jgi$gene, "_unknown_HOLDOUTTEST")
names(jgi_sqs)

# Pull in Marc's training data
nrps <- readAAStringSet("data/sp2.adomainsD.faa")
monomers <- word(names(nrps), 2, sep = "\t")
ids <- gsub("\\.", "_", word(names(nrps), 3, sep = "\t"))
ids
names(nrps) <- paste0(ids, "_", monomers, "_aminoacid_NRPS")
# names(nrps) <- paste0(monomers, "_", names(nrps), "_aminoacid_NRPS")
names(nrps) <- gsub("\\\t", "_", names(nrps))
names(nrps)

mondf <- data.frame(cbind(monomers, names(nrps)), stringsAsFactors = F) %>%
  dplyr::group_by(monomers) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
mondf$monomers

# write.csv(matrix(mondf$monomers, ncol = 13), 'data/monomers_table.csv')
table(mondf$monomers) # take one from each monomer
nrps_rep <- nrps[names(nrps) %in% mondf$V2]

# Search for duplicates with training sqs
cmbnd <- AAStringSet(c(jgi_sqs, training_sqs, nrps_rep))
nodups <- cmbnd[!duplicated(cmbnd)]
length(nodups)
head(names(nodups))
table(word(names(cmbnd), -2, sep = "_"))
length(cmbnd)
names(cmbnd)
writeXStringSet(cmbnd, "data/756_full_ANL_training_set.fasta")

# Align sequences structurally using DECIPHER package 
aa.al <- AlignSeqs(nodups)

# Adjust alignment
aa.adj <- AdjustAlignment(aa.al)
length(aa.adj)
writeXStringSet(aa.adj, "data/725_training_seqs_aligned_adjusted.fasta")

# Trim all columns that have a gap using the trimAl tool
# http://trimal.cgenomics.org/
# ./trimal -in ~/Documents/Wageningen_UR/github/mibig_training_set_build/data/710_training_seqs_aligned_adjusted.fasta -gt 0.9 -out ~/Documents/Wageningen_UR/github/mibig_training_set_build/data/710_trainign_seqs_trimal_90.fasta

# Read in the trimmed alignment
trimaln <- readAAStringSet("data/725_training_seqs_trimal_90.fasta")
width(trimaln) # 398
# BrowseSeqs(trimaln)
# names(trimaln)[grep("cyclosporin", names(trimaln))]
table(word(names(trimaln), -2, sep = "_"))

rdaln <- read.alignment(file = "data/725_training_seqs_trimal_90.fasta", format = "fasta")
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
write.csv(aap, "data/725_seqs_1544_feats_for_supervised.csv", row.names=rownames(aap), quote = F)



