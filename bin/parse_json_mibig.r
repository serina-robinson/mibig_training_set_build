# Install packages
pacman::p_load(
  "tidyverse",
  "data.table",
  "ape",
  "ggtree",
  "ggplot2",
  "pamr",
  "RColorBrewer",
  "phangorn", 
  "ggthemes",
  "usethis",
  "Biostrings",
  "DECIPHER",
  "officer",
  "rvg",
  "drc",
  "stringr",
  "jsonlite",
  "readxl")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build/")

# Look at mibig clustering results
mibig <- readAAStringSet("data/369_MIBiG_non_NRPS_domains.fasta")
head(mibig)

# Trim anything that is longer than 700 aa
# mibig_short <- mibig[width(mibig) < 700]
# mibig_long <- mibig[width(mibig) >= 700]
# names(mibig_long)

# Trim anything that has peptide synthase in the name
no_nrps <- mibig[-grep("peptide", names(mibig))]

# Parse the names
prot_names <- word(names(no_nrps), 6, sep = "\\|")
bgc_names <- word(names(no_nrps), 1, sep = "\\|")
bgc_acc <- word(names(no_nrps), 7, sep = "\\|")
names(no_nrps)
# sort(table(prot_names), decreasing = T)[1:20]
prot_df <- data.frame(cbind(as.character(no_nrps), bgc_acc, prot_names, bgc_names))
colnames(prot_df) <- c("seqs", "acc", "prot_names", "bgcs")
head(prot_df)
# write.csv(cbind(as.character(no_nrps), prot_names), "output/mibig_prot_names.csv")

bgc_ids <- word(names(no_nrps), 1, sep = "\\|")
length(bgc_ids) # 355

# Read in the JSON files for all proteins
allfils <- list.files("~/Documents/Wageningen_UR/github/mibig_training_set_build/data/",full.names = T)
bgcs <- allfils[grep(paste0(bgc_ids, collapse = "|"), allfils)]
bgcs #only 297 had .json files

json_list <- list()
pub_list <- list()
for(i in 1:length(bgcs)) {
  json_list[[i]] <- fromJSON(bgcs[i], simplifyDataFrame = TRUE)
  names(json_list)[i] <- bgcs[i]
  tmp <- json_list[[i]]$general_params$publications
  tmp_spl <- trimws(unlist(str_split(tmp, pattern = ",")))
  pub_list[[i]] <- tmp_spl
}
# str(json_list[[2]])

names(pub_list) <- gsub(".json", "", word(bgcs, 2, sep = "//"))
names(json_list) <- names(pub_list)

to_pull <- unlist(pub_list)
to_pull_un <- unique(to_pull)

to_pull_no_sp_chars <- to_pull_un[-grep(paste0(c("-","unpublished", "10\\.", "\\?"), collapse = "|"), to_pull_un)]
to_pull_final <- to_pull_no_sp_chars[to_pull_no_sp_chars!=""]

finall <- pub_list %>%
  map(~paste0(., collapse = ";"))
finall

unl <- unlist(finall)

# write.table(to_pull_final, quote = F, row.names = F, col.names = F, "output/mibig_pmids_to_pull.txt")

dtf <- data.frame(cbind(unl, names(unl)), stringsAsFactors = F) %>%
  mutate(pmid_last = word(unl, -1, sep = ";"))


colnames(dtf) <- c("pmid", "bgcs", "pmid_last")
head(dtf)
# write_csv(dtf, "output/pmids_for_mibig.csv")


# Read in the downloaded Pubmed data
dat <- read_csv("data/mibig_pubmed_result_20190225.csv") %>%
  janitor::clean_names() %>%
  dplyr::select(title, description, details, entrez_uid) %>%
  dplyr::mutate(pmid_last = as.character(entrez_uid))

# Create a data frame with mibig information
tail(str(json_list))
cmpnd <- list()

chem_struct <- list()
cmpnd <- list()
gene_func <- list()
for(i in 1:length(bgc_join$bgcs)) {
  ind <- grep(bgc_join$bgcs[i], names(json_list))
  cmpnd[[i]] <- json_list[[ind]]$general_params$compounds$compound
  names(cmpnd)[i]<- bgc_join$bgcs[i]
  chem_struct[[i]] <- json_list[[ind]]$general_params$compounds$chem_struct
  gene_func[[i]] <- json_list[[ind]]$general_params$genes$gene$gene_function
  # names(chem_struct)[i] <- bgc_join$bgcs[i] 
}
gene_func

# length(cmpnd)
# names(cmpnd)

cmpnd_df <- data.frame(cbind(names(cmpnd), cmpnd, chem_struct), stringsAsFactors = F)
colnames(cmpnd_df)[1] <- c("bgcs")
cmpnd_df$bgcs <- as.character(unlist(cmpnd_df$bgcs))
class(cmpnd_df$bgcs)
unlist(cmpnd_df$bgcs)

# Merge with the BGC/PMID data
bgc_join <- dat %>%
  inner_join(., dtf, by = "pmid_last") %>%
  inner_join(., cmpnd_df, by = "bgcs") %>%
  unnest(cmpnd) %>%
  left_join(prot_df, ., by = "bgcs") %>%
  distinct(., acc, .keep_all = TRUE) %>%
  write_delim(., "data/mibig_merged_with_pmids.txt", delim = "\t")

writeLines(unique(na.omit(bgc_join$pmid_last)), sep = ",") # to import into Zotero
# pull PDFs when possible

# Read in the existing training data
train_dat <- read_excel("../adenylate_machine_learning/data/anl_training_set_updated_20190215.xlsx")

joint <- train_dat[train_dat$aa_seq %in% bgc_join$seqs,]
head(joint)

# Create final training set
bgc_edited <- read_excel("data/mibig_training_set_manually_edited.xlsx") %>%
    dplyr::filter(!seqs %in% train_dat$aa_seq) %>%
    write_excel_csv(., "data/mibig_training_set_current_20190226.xlsx")

# remove the 5 which are already in the training set

