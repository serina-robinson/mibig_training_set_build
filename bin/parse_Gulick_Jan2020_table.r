# Install packages
pacman::p_load("Biostrings", "bio3d", "readxl", "tidyverse", "janitor")

#Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in Gulick table
gulick <- read_excel("data/20201804_Gulick_ANL_table_parsed.xlsx", col_names = F) %>%
  janitor::clean_names()

pdbr <- list()
for(i in 1:nrow(pdbs)) {
  pdbr[[i]] <- as.character(pdbs[i,1]) 
}


pdbs <- unlist(pdbr)

pdb_pull1 <- na.omit(unlist(str_extract_all(pdbs, "[:digit:][A-Z0-9]{3}")))
pdb_pull2 <- unique(pdb_pull1)
pdb_pull2
pdb_pull3 <- sort(pdb_pull2)

# The last structure included in the original training set was 5VPV
ind <- grep("5VPV", pdb_pull3)
pdb_pull4 <- pdb_pull3[(ind+1):length(pdb_pull3)]
pdb_pull4
writeLines(pdb_pull4, sep = ",")

# PDB ids I want to keep
# # 5X8F,6H1B,6HDY,6IJB,6AAA,6ABH
# rawdat <- read_csv('data/rcsb_pdb_sequence_20200418205518.csv') %>%
#   janitor::clean_names() %>%
#   dplyr::filter(grepl("5X8F|6H1B|6HDY|6IJB|6AAA|6ABH", pdb_id))

# 6H1B is only 33%
# 6HDY is 33.99%
# 6IJB is 81% to Ruegeria pomeroyi C2 through C5
# 6AAA is 86% to Lampyris luciferase
# 6ABH and 5X8F are removed

# rawdat$sequence[3] <- gsub("MASHHHHHHSG", "", rawdat$sequence[3])
# towrite <- AAStringSet(rawdat$sequence)
# names(towrite) <- rawdat$pdb_id
# writeXStringSet(towrite, "output/6_new_PDB_struct_seqs.fasta")

# To test AdenylPred
rawdat <- read_csv('data/rcsb_pdb_sequence_20200418205518.csv') %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl("6H1B|6HDY|6IJB|6AAA", pdb_id))
rawdat$sequence[2] <- gsub("MASHHHHHHSG", "", rawdat$sequence[2])
towrite <- AAStringSet(rawdat$sequence)
names(towrite) <- rawdat$pdb_id
length(towrite)
names(towrite)

writeXStringSet(towrite, "output/4_new_PDB_struct_seqs.fasta")
write_csv(rawdat, "output/4_new_PDB_struct_features.csv")
                

# Read in the table
tab <- read_excel("output/20200415_JBC_MIBiG_validation_set.xlsx") #%>%
tab$aa_seq
sqs <- readAAStringSet("data/MIBiG_accessions_since_2019.fasta")
pdbs <- readAAStringSet("output/4_new_PDB_struct_seqs.fasta")
pdbs
rue <- readAAStringSet("data/Ruegeria_pomeroyi.fasta")
names(pdbs)

mibig <- sqs[grep(paste0(tab$acc, collapse = "|"), names(sqs))]
tab$aa_seq
total <- AAStringSet(c(AAStringSet(tab$aa_seq[!is.na(tab$aa_seq)]), mibig, rue, pdbs))
names(total)[1:length(tab$aa_seq[!is.na(tab$aa_seq)])] <- paste0(tab$query_name[!is.na(tab$aa_seq)],
                                                                 "_",
                                                                 tab$organism[!is.na(tab$aa_seq)],
                                                                 "_",
                                                                 tab$true_or_highly_likely_substrates[!is.na(tab$aa_seq)])

dedup <- total[!duplicated(total)]
names(dedup) <- gsub("\\||\\;|\\?|\\-|\\[|\\]|\\+| |\\,|\\(|\\)", "_", names(dedup))
writeXStringSet(dedup, "output/30_holdout_sequences_for_validation.fasta")

# Read in the 1553 training seqs
allsqs <- readAAStringSet("data/1761_combined_training_full_length_sqs.faa")
