# Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Get the tuning parameters
nnet_ml <- readRDS("data/model_comparisons/nnet_NRPS_broad_groups_repeatedCV_20190306.rds")
nnet_ml$bestTune

# Read in the 34aa HMM aligned sequences
aa34 <- readAAStringSet("~/Documents/Wageningen_UR/github/check_a_dom_specificity/data/sp2_adoms_34extract_hmmalign.faa")
name_key <- readAAStringSet("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/data/sp2_34extract_names_fixed_large_grps.faa")
names(aa34) <- names(name_key)
table(word(names(aa34), sep = "_", -1))

# Remove all duplicates
dedup <- aa34[!duplicated(aa34)]
length(dedup) # 843 
subst <- word(names(dedup), sep = "_", 2)
subst <- gsub("\\|", "\\.", subst)
subst <- gsub("\\-", "\\.", subst)
names(dedup) <- gsub("\\|", "\\.", names(dedup))
names(dedup) <- gsub("\\-", "\\.", names(dedup))
table(subst)
length(subst)
# writeXStringSet(dedup, "data/843_34extract_hmmalign.faa")

dtf <- data.frame(names(dedup), dedup, subst, stringsAsFactors = F) %>%
  dplyr::filter(!grepl("reject", names.dedup.)) %>%
  dplyr::filter(!grepl("-", dedup)) %>%
  group_by(subst) %>%
  add_count(subst) %>%
  dplyr::filter(n > 8) 

dim(dtf) # 652
table(dtf$subst)

# Remove all sequences with specificity less than 8
dedup2 <- AAStringSet(dtf$dedup)
names(dedup2) <- dtf$names.dedup.
writeXStringSet(dedup2, "data/652_monomers_deduplicated.faa")

# First do large substrate groups
rdaln <- read.alignment("data/652_monomers_deduplicated.faa", format = "fasta")
str(rdaln)

# Convert to a vector of physicochemical properties
source("src/convert_aln_15aap.r")
aa <- convert_aln_15aap(rdaln) #15 physicochemical properties
aadf <- data.frame(aa,stringsAsFactors = F)
dim(aadf) # 843 x 510
rownames(aadf)
# write.csv(aadf, paste0("data/652_aa_seqs_510_feats_for_supervised_monomers_20190406.csv"), row.names=rownames(aadf), quote = F)

# Read in the data
rawdat <- read_csv("data/652_aa_seqs_510_feats_for_supervised_monomers_20190406.csv", skip_empty_rows = T)

rawdat <- data.frame(cbind(rawdat$X1), scale(rawdat[,2:ncol(rawdat)]), stringsAsFactors = F)
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Remove invariable columns


# Remove the holdout test predictions
# dat <- rawdat[-grep("reject", rawdat$nms),] # 658 observations
# dim(dat)

# Set seed 
set.seed(20190506)
dim(rawdat)
# dat <- rawdat[!duplicated(rawdat),]
dat <- rawdat %>%
  dplyr::select(-RADA880108_27)
dim(dat)
table(dat$clf)


# Split into test and training data
dat_split <- initial_split(dat, strata = "clf")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
y_train <- as.factor(dat_train$clf)
y_test <- as.factor(dat_test$clf)
table(y_test)


# Make a data frame for prediction
df_train <- data.frame(x_train, stringsAsFactors = F, row.names = dat_train$nms)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)



nnet_grid <- expand.grid(.decay = c(0.5, 1e-3, 1e-7),
                         .size = c(3, 5, 10, 20))

nnet_grid <- expand.grid(.decay = c(0.5),
                         .size = c(10))

nnet_mod <- train(
  x = df_train,
  y = y_train,
  method = "nnet",
  MaxNWts = 100000,
  # size = 1,
  tuneGrid = nnet_grid,
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nnet_mod)
nnet_mod$results

saveRDS(nnet_mod, "data/model_comparisons/nnet_NRPS_monomers_quick_20190406.rds")
nnet_ml <- readRDS("data/model_comparisons/nnet_NRPS_monomers_quick_20190406.rds")
nnet_pred <- predict(nnet_ml, newdata = form_test)
cm_nnet <- confusionMatrix(nnet_pred, as.factor(dat_test$clf))
cm_nnet
0.6933 - 0.6163
