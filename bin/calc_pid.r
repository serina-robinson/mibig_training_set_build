# Install packages
pacman::p_load("DECIPHER", "tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the 575 training set seqs
training_sqs <- readAAStringSet('data/575_training_seqs_including_mibig.fasta')
table(word(names(training_sqs), -2, sep = "_"))
head(training_sqs)

# Read in the JGI seq
jgi <- jgi <- read_excel("data/Proposal_503776_DOE_JGI_genes_final.xlsx") %>%
  janitor::clean_names() %>%
  dplyr::filter(grepl("OleC", role)) %>%
  dplyr::filter(grepl("Bacter", gene))
bact <- AAStringSet(jgi$seq)
names(bact) <- jgi$gene

# Calculate pairwise identities
pairws <- list()
for(i in 1:length(training_sqs)) {
  tmp <- pairwiseAlignment(bact, training_sqs[i])
  pairws[[i]] <- pid(tmp, type="PID1") 
}
head(pairws)
names(pairws) <- names(training_sqs)

# Find the maximum pid
max_pid <- sort(unlist(pairws), decreasing = T)
max_pid[1:20]

# What is the substrate of the maximum



