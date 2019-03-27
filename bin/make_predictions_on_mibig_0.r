## Install packages
pacman::p_load("DECIPHER", "readxl", "data.table", "Biostrings", "caret", "bgafun", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the updated mibig
mib <- read_excel("data/mibig_training_set_manually_edited_20192503.xlsx") %>%
  dplyr::filter(confidence < 1) %>%
  mutate(substrate_group_tr = str_replace_all(substrate_group, "_", "")) %>%
  mutate(sqnams_tr = paste0(bgcs, "_", word(cmpnd, 1, sep = "_"), "_", acc, "_", substrate_group_tr, "_", functional_class))


# Extract the 34 aa signatures
toextract <- AAStringSet(mib$aa_seq)
names(toextract) <-mib$sqnams_tr
writeXStringSet(toextract, '../sandpuma2_serina/flat/185_mibig_low_confidence.fa')

# Extract 34 aa signatures
sigs <- fread('../sandpuma2_serina/flat/185_low_confidence_34_aa_sigs.fa', data.table = F, header = F)
sigaa <- AAStringSet(sigs[,1])
names(sigaa) <- mib$sqnams_tr
writeXStringSet(sigaa, 'data/185_mibig_34_aa_sigs.fa')

# Convert 34 aa signatures to 
rdaln <- read.alignment(file = 'data/185_mibig_34_aa_sigs.fa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)
aa <- bgafun::convert_aln_AAP(rdaln) #5 physicochemical properties

aadf <- data.frame(aa, stringsAsFactors = F)

aap <- aadf %>%
  dplyr::mutate(nms = rownames(.)) %>%
  dplyr::select(-contains("D")) %>%
  dplyr::filter(!grepl("ERROR", X1A))

rownames(aap) <- aap$nms
aap <- aap %>%
  dplyr::select(-nms)
colnames(aap) <- gsub("^X","",colnames(aap))

colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap))
numfeats <- length(colnames(aap)) # 1096
aap[,1]

# Read in the ml model
mod <- readRDS("output/rf_small_subst_grp_20192603.rds")
makpred <- predict(mod, newdata = aap)
makpred

# Write the predictions to file
preds <- write_csv(data.frame(cbind(rownames(aap), as.character(makpred))), "output/185_mibig_and_preds.csv")
