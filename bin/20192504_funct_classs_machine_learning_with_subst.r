## Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190403)

# Read in the data
rawdat <- read_csv("data/703_training_sqs_with_loop_extracted.csv")

table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
head(rawdat)
# rawdat$nms
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Try downsampling the ARYL and NRPS sequecnes to about 80 
# 4-coumarate sequences in the aryl class
sample(which(rawdat$clf == "NRPS"), 80)
sample(which(rawdat$clf == "ARYL"), 80)

rawdat_tr <- rawdat %>%
  dplyr::filter(stringr::str_detect(pattern = "ARYL|NRPS", clf, negate = T)) %>%
  bind_rows(., rawdat[sample(which(rawdat$clf == "NRPS"), 80),]) %>%
  bind_rows(., rawdat[sample(which(rawdat$clf == "ARYL"), 80),])
dim(rawdat_tr)

# Remove the holdout test predictions
dat <- rawdat_tr[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat_tr$nms),] # 658 observations

dat_split <- initial_split(dat, prop = 3/4, strata = "clf")
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
# tunegrid <- expand.grid(.mtry = c(10, as.integer(sqrt(ncol(x_train))), 140), .min.node.size = 1)

# tunegrid

rf <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
 #  tuneGrid = tunegrid,
  verbose = TRUE,
  # max.depth = 20,
  importance = "permutation")
getTrainPerf(rf)
rf$results[which.max(rf$results$Accuracy),]

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
cm_rf

sink("data/model_comparisons/cm_rf_functional_class_20192604.txt")
cm_rf
sink()

dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_Precision_Recall_Accuracy_Functional_class_2012604.csv", row.names = T)

# Try with weights
classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
classwts
dtf <- data.frame(cbind(table(dat$clf), classwts))
colnames(dtf) <- c("class_size", "class_weights")
# write.csv(dtf, "output/small_sub_grp_class_sizes_and_weights.csv", quote = F, row.names = T)

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

rf_wt_pred <- predict(rf_weighted, newdata = form_test)
cm_rf <- confusionMatrix(rf_wt_pred, as.factor(dat_test$clf))
cm_rf

sink("data/model_comparisons/cm_rf_wt_pred_20192604.txt")
cm_rf
sink()

dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_Precision_Recall_Accuracy_Functional_class_2012604_weighted.csv", row.names = T)


### RF weighted with a max depth of 15
rf_15 <- caret::train(
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
  max.depth = 15,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted)

rf_15_pred <- predict(rf_15, newdata = form_test)
cm_rf <- confusionMatrix(rf_15_pred, as.factor(dat_test$clf))
cm_df <- data.frame(melt(cm_rf$table), stringsAsFactors = F)
head(cm_df)

pdf("data/cm_rf_funct_class_confusion_matrix.pdf")
ggplot(cm_df, aes(Reference, Prediction)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()


sink("data/model_comparisons/cm_rf_15_pred_20192604.txt")
cm_rf
sink()

dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select
write.csv(dtl_feat_select, "data/model_comparisons/rf_15_AA34_Precision_Recall_Accuracy_Functional_class_2012604_weighted.csv", row.names = T)

