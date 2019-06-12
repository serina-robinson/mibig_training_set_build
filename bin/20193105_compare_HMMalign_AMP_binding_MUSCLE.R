# Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190304)

# Read in the two alignments 
musc <- read.alignment("data/20192305_MUSCLE_34aa_extracted.faa", format = "fasta")
hmm_adom <- read.alignment("data/20192305_1693_small_substrate_grp_for_rf_34extracted_a_dom_hmm.faa", format = "fasta")
hmm_amp <- read.alignment("data/20192305_1693_small_substrate_grp_34extracted_AMP_binding.fasta", format = "fasta")

# Convert to a feature matrix
# source("src/convert_aln_15aap.r")
# musc_df <- convert_aln_15aap(musc)
# rownames(musc_df) <- musc$nam
# write.csv(data.frame(musc_df, stringsAsFactors = F), "output/20192305_MUSCLE_15aa_properties.csv", quote = F, row.names = T)
# 
# # Convert HMMAlign adom results to feature matrix
# hmm_adom_df <- convert_aln_15aap(hmm_adom)
# rownames(hmm_adom_df) <- hmm_adom$nam
# write.csv(data.frame(hmm_adom_df, stringsAsFactors = F), "output/20192305_HMMAlign_15aa_a_dom_properties.csv", quote = F, row.names = T)
# 
# # Convert HMMAlign AMP-binding results to feature matrix
# hmm_amp_df <- convert_aln_15aap(hmm_amp)
# dim(hmm_amp_df)
# rownames(hmm_amp_df) <- hmm_amp$nam
# write.csv(data.frame(hmm_amp_df, stringsAsFactors = F), "output/20192305_HMMAlign_15aa_amp_binding_properties.csv", quote = F, row.names = T)

musc_df <- read_csv("output/20192305_MUSCLE_15aa_properties.csv")
hmm_adom_df <- read_csv("output/20192305_HMMAlign_15aa_a_dom_properties.csv")
hmm_amp_df <- read_csv("output/20192305_HMMAlign_15aa_amp_binding_properties.csv")
df_list <- list(musc_df, hmm_adom_df, hmm_amp_df)
length(df_list)

rf_results_list <- list()
for(i in 1:length(df_list)) {
  
  tmp <- df_list[[i]]
  rawdat <- data.frame(tmp[!duplicated(tmp[,2:ncol(tmp)]),], stringsAsFactors = F)
  
  table(duplicated(rawdat[,2:ncol(rawdat)])) # check none are duplicated
  colnames(rawdat)[1] <- "nms"
  rawdat$clf <- word(rawdat$nms, 2, sep = "_")
  table(rawdat$clf)
  
  # Remove rejects
  dat  <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "reject"), collapse = "|"), rawdat$nms),]

  # Train a random forest model with optimal parameters, no max depth
  x_train <- dat[,!colnames(dat) %in% c("nms", "clf")]
  y_train <- dat$clf
  
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$nms)
  
  # Compare the different models
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
  
    rf_results_list[[i]] <- rf_full_ss
}

# Compare for small substrate groups

# MUSCLE alignment is pretty good
rf_results_list[[1]]$results
# mtry splitrule min.node.size  Accuracy     Kappa AccuracySD    KappaSD
# 1   22      gini             1 0.8131467 0.7944241 0.02886966 0.03181204

# A-dom alignment is on par
rf_results_list[[2]]$results
# mtry splitrule min.node.size  Accuracy     Kappa AccuracySD    KappaSD
# 1   22      gini             1 0.8188931 0.8007196 0.03024923 0.03325956

# AMP-Binding HMM
rf_results_list[[3]]$results # Significantly worse performance than the other two
# mtry splitrule min.node.size  Accuracy     Kappa AccuracySD    KappaSD
# 1   22      gini             1 0.7679378 0.7426659 0.03342423 0.03711207


rf_results_list <- list()
for(i in 1:length(df_list)) {
  
  tmp <- df_list[[i]]
  rawdat <- data.frame(tmp[!duplicated(tmp[,2:ncol(tmp)]),], stringsAsFactors = F)
  
  table(duplicated(rawdat[,2:ncol(rawdat)])) # check none are duplicated
  colnames(rawdat)[1] <- "nms"
  rawdat$clf <- word(rawdat$nms, 2, sep = "_")
  table(rawdat$clf)
  
  # Remove rejects
  dat  <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "reject"), collapse = "|"), rawdat$nms),]
  
  # Train a random forest model with optimal parameters, no max depth
  x_train <- dat[,!colnames(dat) %in% c("nms", "clf")]
  y_train <- dat$clf
  
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$nms)
  
  # Compare the different models
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
  
  rf_results_list[[i]] <- rf_full_ss
}

saveRDS(rf_results_list[[1]], "data/model_comparisons/rf_MUSCLE_aligned_ss_fulltrain_20193105.rds")
saveRDS(rf_results_list[[2]], "data/model_comparisons/rf_A_dom_aligned_ss_fulltrain_20193105.rds")
saveRDS(rf_results_list[[3]], "data/model_comparisons/rf_AMP_binding_aligned_ss_fulltrain_20193105.rds")


rf_results_fc <- list()
for(i in 1:length(df_list)) {
  
  tmp <- df_list[[i]]
  rawdat <- data.frame(tmp[!duplicated(tmp[,2:ncol(tmp)]),], stringsAsFactors = F)
  
  table(duplicated(rawdat[,2:ncol(rawdat)])) # check none are duplicated
  colnames(rawdat)[1] <- "nms"
  rawdat$clf <- word(rawdat$nms, 3, sep = "_")
  table(rawdat$clf)
  
  # Remove rejects
  dat  <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "reject"), collapse = "|"), rawdat$nms),]
  
  # Train a random forest model with optimal parameters, no max depth
  x_train <- dat[,!colnames(dat) %in% c("nms", "clf")]
  y_train <- dat$clf
  
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat$nms)
  
  # Compare the different models
  tunegrid <- expand.grid(.mtry = as.integer(sqrt(ncol(x_train))), .splitrule = 'gini', .min.node.size = 1)
  
  rf_full_fc <- caret::train(
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
  
  rf_results_fc[[i]] <- rf_full_fc
}
rf_results_fc[[1]]$results
rf_results_fc[[2]]$results
rf_results_fc[[3]]$results

saveRDS(rf_results_fc[[1]], "data/model_comparisons/rf_MUSCLE_aligned_fc_fulltrain_20193105.rds")
saveRDS(rf_results_fc[[2]], "data/model_comparisons/rf_A_dom_aligned_fc_fulltrain_20193105.rds")
saveRDS(rf_results_fc[[3]], "data/model_comparisons/rf_AMP_binding_aligned_fc_fulltrain_20193105.rds")

