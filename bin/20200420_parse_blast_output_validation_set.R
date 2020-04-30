# Install packages
pacman::p_load("Biostrings", "bio3d", "readxl", "tidyverse", 
               "janitor", "dplyr", "data.table")

#Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the table
tab <- read_excel("output/20200421_JBC_MIBiG_validation_set_final.xlsx") %>%
  dplyr::slice(1:41)
tab$acc

tabr <- AAStringSet(tab$aa_seq[tab$aa_seq != "NA"])
tabr
names(tabr) <- tab$acc[tab$aa_seq != "NA"]

sqs <- readAAStringSet("output/full_length_MIBiG_new_BGCs_Adomains_not_extracted.fasta")
tab$acc

mibig <- sqs[grep(paste0(tab$acc, collapse = "|"), names(sqs))]
names(mibig)
names(mibig) <- word(names(mibig), sep = "\\|", -1)


pdbs <- readAAStringSet("output/4_new_PDB_struct_seqs.fasta")
names(pdbs) <- paste0(names(pdbs), "_A")

rue <- readAAStringSet("data/Ruegeria_pomeroyi.fasta")
names(rue) <- "WP_011047771.1"

slividans <- readAAStringSet("data/20200420_S_lividans_seqs.fasta")
names(slividans) <- word(names(slividans), sep = " ", 1)

total <- AAStringSet(c(tabr, mibig, rue, pdbs, slividans))
total

dedup <- total[!duplicated(total)]



#writeXStringSet(dedup, "output/40_holdout_sequences_for_validation_names_fixed.fasta")

# Read in the BLAST hits
# 
blast <- fread("blast/41_holdout_max_blast_hits.tsv", data.table = F, header = F) %>%
  dplyr::filter(!grepl("HOLDOUTTEST|OTHER", V2)) %>%
  dplyr::group_by(V1) %>%
  dplyr::select(V1, V2, V3) %>%
  slice(1)
blast$V3

colnames(blast) <- c("query", "hit", "perc_id")
dim(blast)

# write_csv(blast, "output/29_hits_with_pid_scores.csv")


# Combine with the training set
jbc <- read_excel('output/20200421_JBC_MIBiG_validation_set_final.xlsx') %>%
  janitor::clean_names() %>%
  dplyr::slice(1:41)

jbc$acc[jbc$acc == "D581_RS0125480"] <- "WP_018960015.1"
jbc$acc[jbc$acc == "VIBHAR_01399"] <- "ABU70375.1"
jbc$acc

setdiff(jbc$acc, blast$query)
setdiff(blast$query, jbc$acc)

comb <- jbc %>%
  dplyr::right_join(blast, ., by = c("query" = "acc"))
colnames(comb)
combfindf <- comb %>%
  dplyr::filter(!is.na(perc_id)) %>%
  dplyr::filter(perc_id != 100.0) #%>%
 # dplyr::select(query, hit, perc_id, true_or_highly_likely_substrates, fc_correct, )

dim(combfindf)
head(combfindf)

write_csv(combfindf, "output/40_validation_seqs_for_labmeeting_20200420.csv")
