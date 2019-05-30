# Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190304)

# Read in the two alignments 
musc <- read.alignment("data/20192305_MUSCLE_34aa_extracted.faa", format = "fasta")
hmm <- read.alignment("data/20192305_1693_small_substrate_grp_for_rf_34extracted_a_dom_hmm.faa", format = "fasta")


# Convert to a feature matrix
source("src/convert_aln_15aap.r")
musc_df <- convert_aln_15aap(musc)
rownames(musc_df) <- musc$nam
write.csv(data.frame(musc_df, stringsAsFactors = F), "output/20192305_MUSCLE_15aa_properties.csv", quote = F, row.names = T)

hmm_df <- convert_aln_15aap(hmm)
rownames(hmm_df)
write.csv(data.frame(hmm_df, stringsAsFactors = F), "output/20192305_HMMAlign_15aa_properties.csv", quote = F, row.names = T)

musc_df <- read_csv("output/20192305_MUSCLE_15aa_properties.csv")


rawdat <- data.frame(musc_df[!duplicated(musc_df[,3:ncol(musc_df)]),], stringsAsFactors = F) %>%
  dplyr::select(-X1_1)
colnames(rawdat)
dim(rawdat) # 1305
table(duplicated(rawdat[,2:ncol(rawdat)])) # none are duplicated
colnames(rawdat)[1] <- "nms"
# colnames(rawdat)[1]
class(rawdat)
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Remove rejects
# dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "amino.acid", "reject"), collapse = "|"), rawdat$nms),]
# dim(dat)
dat <- rawdat
colnames(dat)

# Train a random forest model with optimal parameters, no max depth
x_train <- dat[,!colnames(dat) %in% c("nms", "clf")]
y_train <- dat$clf

form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$nms)


tunegrid <- expand.grid(.mtry = as.integer(sqrt(ncol(x_train))), .splitrule = 'gini', .min.node.size = 1)

rf_full_ss <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  tuneGrid = tunegrid,
  verbose = TRUE,
  importance = "permutation")
rf_full_ss$results
