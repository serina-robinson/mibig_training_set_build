## Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20191304)

# Read in the training set
grp_comb <- readAAStringSet("output/20191604_small_substrate_grp_for_rf_duplicates_included_py.fasta")
summary(width(rawdat))

# Convert to seq features
source("src/convert_seq_15aap.r")
res <- strsplit(as.character(grp_comb), "")
nlist <- lapply(1:length(grp_comb), function(x) convert_seq_15aap(as.character(res[[x]])))
extract_34_df <- data.frame(matrix(unlist(nlist), nrow = length(res), byrow=T), stringsAsFactors=FALSE)

colnames(extract_34_df) <- as.character(fread("data/feature_names.txt", data.table = F)[,1])
dim(extract_34_df) # 2344 predictions by 510
rownames(extract_34_df) <- names(grp_comb)
head(extract_34_df)
write.csv(extract_34_df, "data/1697_seqs_510_features_centered_scaled.csv", quote = F, row.names = T)

# Set seed 
set.seed(20191304)

# Read in the data
rawdat <- read_csv("data/1697_seqs_510_features_centered_scaled.csv")
colnames(rawdat)[1] <- "nms"
head(rawdat)
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid"), collapse = "|"), rawdat$nms),] # 658 observations
table(dat$clf)

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

# Try without weights
tunegrid <- expand.grid(.splitrule = "gini", .mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)

tunegrid

rf <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  tuneGrid = tunegrid,
  verbose = TRUE,
  max.depth = 9,
  importance = "permutation")
getTrainPerf(rf)


# Try with weights
classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))

dtf <- data.frame(cbind(table(dat$clf), classwts))
colnames(dtf) <- c("class_size", "class_weights")
write.csv(dtf, "output/small_sub_grp_class_sizes_and_weights.csv", quote = F, row.names = T)

rf_weighted <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 500,
  tuneGrid = tunegrid,
  verbose = TRUE,
  max.depth = 9,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted)

rf_15 <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  tuneGrid = tunegrid,
  verbose = TRUE,
  max.depth = 15,
  importance = "permutation")
getTrainPerf(rf_15)
getTrainPerf(rf)

rf_12 <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  tuneGrid = tunegrid,
  verbose = TRUE,
  max.depth = 12,
  importance = "permutation")
getTrainPerf(rf_12)
getTrainPerf(rf_15)

# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}
saveRDS(rf_15, "data/model_comparisons/rf_15_xvalidated_small_subst_grp_201901604.rds")

# ROC curve is not looking good...
approx_roc_curve(rf_15, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# Try prediction
# rf_ml <- ranger(y_train ~., data = form_train, num.trees = 500, splitrule = as.character(rf$bestTune$splitrule), 
#                 mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
#                 importance = "permutation")
# saveRDS(rf_ml, "data/model_comparisons/rf_small_subst_group_20191604.rds")
# rf_ml <- readRDS("data/model_comparisons/rf_small_subst_group_20190204.rds")

rf_pred <- predict(rf_15, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))

sink("data/model_comparisons/cm_rf_15aa_functional_class_20191604.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T), 
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp

pdf("data/model_comparisons/rf_feat_select_var_imp_small_subst_group_20191604.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()


dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_Precision_Recall_Accuracy_Substrate_group_20191604.csv", row.names = T)

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
ggsave("data/model_comparisons/aa34_rf_conf_matrices_subst_grp_20191604.jpeg", height=7, width=7, units='in')


