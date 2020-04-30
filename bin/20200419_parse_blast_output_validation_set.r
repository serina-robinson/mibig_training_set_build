# Install packages
pacman::p_load("Biostrings", "bio3d", "readxl", "tidyverse", "janitor", "dplyr")

#Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")


# Read in the table
# tab <- read_excel("output/20200415_JBC_MIBiG_validation_set.xlsx") 
# tabr <- AAStringSet(tab$aa_seq[!is.na(tab$aa_seq)])
# names(tabr) <- tab$acc[!is.na(tab$aa_seq)]
# sqs <- readAAStringSet("data/MIBiG_accessions_since_2019.fasta")
# mibig <- sqs[grep(paste0(tab$acc, collapse = "|"), names(sqs))]
# names(mibig) <- word(names(mibig), sep = "\\|", -1)
# 
# pdbs <- readAAStringSet("output/4_new_PDB_struct_seqs.fasta")
# names(pdbs) <- paste0(names(pdbs), "_A")
# 
# rue <- readAAStringSet("data/Ruegeria_pomeroyi.fasta")
# names(rue) <- "WP_011047771.1"
# 
# total <- AAStringSet(c(tabr, mibig, rue, pdbs))
# 
# dedup <- total[!duplicated(total)]
# names(dedup)
# writeXStringSet(dedup, "output/30_holdout_sequences_for_validation_names_fixed.fasta")

# Read in the BLAST hits
# RID-9U1PSZ15114
blast <- read_csv("data/9U1PSZ15114-Alignment-HitTable-namsfixed.csv", col_names = F) %>%
  dplyr::filter(!grepl("HOLDOUTTEST|OTHER", X2)) %>%
  dplyr::group_by(X1) %>%
  dplyr::select(X1, X2, X3) %>%
  slice(1)
colnames(blast) <- c("query", "hit", "perc_id")
blast$query

# write_csv(blast, "output/29_hits_with_pid_scores.csv")
setdiff(jbc$acc, blast$query)
setdiff(blast$query, jbc$acc)

# Combine with the training set
jbc <- read_excel('output/20200415_JBC_MIBiG_validation_set.xlsx') %>%
  janitor::clean_names() 

jbc$acc[jbc$acc == "D581_RS0125480"] <- "WP_018960015.1"
jbc$acc[jbc$acc == "VIBHAR_01399"] <- "ABU70375.1"

comb <- jbc %>%
  dplyr::right_join(blast, ., by = c("query" = "acc"))
colnames(comb)
combfindf <- comb %>%
  dplyr::filter(!is.na(perc_id))
dim(combfindf)
write_csv(combfindf, "output/29_validation_seqs_for_labmeeting_20200419.csv")
