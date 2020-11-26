# Install packages
pacman::p_load("tidyverse", "readxl", "Biostrings", "data.table")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the dataset
dat <- read_csv("data/F4SZZFCG01N-Alignment-HitTable-11296_sequences.csv", col_names = F)
summary(dat$X3) 
head(dat)
length(unique(dat$X2))
nrow(dat)
accs <- dat %>% 
  dplyr::pull(2) %>%
  unique()
length(accs)
write_csv(data.frame(accs), "output/10000_accs_for_batch_entrez", col_names = F)

dat <- readAAStringSet("data/10000_BLS.fasta")
dat_dedup <- dat[!duplicated(dat)]
dat_deduplength(dat)


datdf <- data.frame(nams = names(dat), seqs = dat)

datdf$nams
gendf <- datdf %>%
  dplyr::mutate(genus = gsub("\\]", "", stringr::word(stringr::word(nams, sep = "\\[", -1), sep = " ", 1))) %>%
  group_by(genus) %>%
  dplyr::slice(1)
gendf$genus
# genaa <- AAStringSet(gendf$seqs[750:849])
# names(genaa) <- gendf$nams[750:849]
# names(genaa)
writeXStringSet(tow, "output/last_100_hits.fasta")

genun <- unique(gendf$genus)

length(genun)
genun

nrow(dat)
tow <- AAStringSet(dat[8824:8923])
names(tow)

# Taxonomy
taxx <- data.table::fread("data/parsed_tax.tsv", sep = "\t", data.table = F) %>%
  janitor::clean_names() %>%
  dplyr::select(-species, -genome_index) %>%
  distinct()

mergdf <- gendf %>%
  dplyr::left_join(., taxx, by = "genus")
length(table(mergdf$phylum)) # 18 different phyla
length(table(mergdf$family)) #191
table(mergdf$family)
# Merge with dataset