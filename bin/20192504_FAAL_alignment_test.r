# Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "readxl", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
comb <- read_excel("data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
faals <- comb %>%
  dplyr::filter(functional_class == "FAAL") %>%
  mutate(sqnams = paste0(acc, "_", organism, "_", small_substrate_group, "_", functional_class))

faal_sqs <- AAStringSet(faals$aa_seq)
names(faal_sqs) <- faals$sqnams
# writeXStringSet(faal_sqs, "data/all_faals.fasta")
# writeXStringSet(faal_sqs[grep("FAA32_MYCTU", names(faal_sqs))], "data/myctu_faal32.fasta")

# Align the FAALs
faal_aln <- AlignSeqs(faal_sqs)
BrowseSeqs(faal_aln)

# Align FAALs with Adomains muscle
ref <- readAAStringSet("data/A_domains_muscle.fasta")

try_al <- AlignProfiles(faal_aln, ref)
BrowseSeqs(try_al)


sqs <- readAAStringSet(input$file1$datapath)
query_fils <- sapply(1:length(sqs), function(x) {tempfile(pattern = "", fileext = ".fasta")})
sapply(1:length(sqs), function(x) {writeXStringSet(sqs[x], query_fils[x])})
extract_34_list <- lapply(1:length(sqs), function(x) { extract_34_aa(query_fils[x]) })
extract_34_df <- data.frame(matrix(unlist(extract_34_list), nrow = length(sqs), byrow=T), 
                            stringsAsFactors=FALSE)
colnames(extract_34_df) <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
