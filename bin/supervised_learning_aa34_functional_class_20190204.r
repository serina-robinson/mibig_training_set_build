## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
rawdat <- readAAStringSet("data/669_training_set_34aa_extracted_20192703.faa")
length(rawdat)
clf <- word(names(rawdat), -1, sep = "_")
subst1<- word(names(rawdat), -2, sep = "_")
subst2<- word(names(rawdat), -3, sep = "_")

# Read in the NRPS training set
aa_dat <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
ran_nums <- sample(1:length(aa_dat), 150, replace = F)
aa_dat <- aa_dat[ran_nums]

names(aa_dat) <- paste0(names(aa_dat), "_aminoacid_NRPS")
clf <- word(names(aa_dat), -1, sep = "_")
subst1<- word(names(aa_dat), -2, sep = "_")
subst2<- word(names(aa_dat), -3, sep = "_")
table(subst2)
aa_rem <- aa_dat[subst2 != "reject"]

# Combine the two
comb <- AAStringSet(c(rawdat, aa_rem))
length(comb)
writeXStringSet(comb, "data/818_aa34_signatures_20190104.fa")

# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/818_aa34_signatures_20190104.fa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)

source("src/convert_aln_15aap.r")
aa <- convert_aln_15aap(rdaln) #15 physicochemical properties
aadf <- data.frame(aa, stringsAsFactors = F)

# Write to CSV file
# write.csv(aadf, paste0("data/818_seqs_510_feats_for_supervised_20190204.csv"), row.names=rownames(aadf), quote = F)

# Read in the data
rawdat <- read_csv("data/818_seqs_510_feats_for_supervised_20190204.csv")
rawdat <- data.frame(cbind(rawdat$X1), scale(rawdat[,2:ncol(rawdat)]), stringsAsFactors = F)
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid"), collapse = "|"), rawdat$nms),] # 658 observations

# Set seed 
set.seed(20190304)

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

####### Random Forest with feature selection ########
# Random Forest using the ranger package for prediction
# Fitting mtry = 2, splitrule = gini, min.node.size = 1 on full training set
rf <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation",
  preProcess = c("center", "scale"))

saveRDS(rf, "data/model_comparisons/rf_xvalidated_functional_class_20190204.rds")

# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}

# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}
saveRDS(rf, "data/model_comparisons/rf_xvalidated_functional_class_20190204.rds")

# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# Try prediction
rf_ml <- ranger(y_train ~., data = form_train, num.trees = 500, splitrule = as.character(rf$bestTune$splitrule), 
                mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")
saveRDS(rf_ml, "data/model_comparisons/rf_functional_class_20190204.rds")
rf_ml <- readRDS("data/model_comparisons/rf_functional_class_20190204.rds")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
cm_rf
sink("data/model_comparisons/cm_rf_15aa_functional_class_20190204.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T), 
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp

pdf("data/model_comparisons/rf_feat_select_var_imp_functional_class_20190204.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()


dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/rf_functional_class_AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_rf$table)
names(cm_list) <- c("rf_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/model_comparisons/aa34_rf_functional_class_conf_matrices_subst_grp_20190204.jpeg", height=7, width=7, units='in')
