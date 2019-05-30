# Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190304)

# Read in the data
rawdat <- read_csv("data/1553_training_sqs_with_loop_extracted.csv")
table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
head(rawdat)
rawdat$nms
rawdat$clf <- word(rawdat$nms, -1, sep = "_")
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat$clf),] # 658 observations
dat_split <- initial_split(dat, prop = 3/4, strata = "clf")
dat_train <- training(dat_split)
dat_test  <- testing(dat_split)
nrow(dat_train)/nrow(dat) # 75 %
table(word(dat$nms, sep = "_", 2))

# Define our response
x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
y_train <- as.factor(dat_train$clf)
y_test <- as.factor(dat_test$clf)
table(y_test)

# Complete dataset for training and testing
form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)

classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
classwts
dtf <- data.frame(cbind(table(dat$clf), classwts))
colnames(dtf) <- c("class_size", "class_weights")

# writev(dtf, "output/small_sub_grp_class_sizes_and_weights.csv", quote = F, row.names = T)

### Try it with max depth 15 and tune grid determined
tunegrid <- expand.grid(.splitrule = "gini", .mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)

rf_weighted <- caret::train(
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
  max.depth = 15,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted)

rf_15_wtd <- rf_weighted
saveRDS(rf_15_wtd, "data/20190305_rf_models/rf_15_wtd_small_sub_grp.rds")
rf_wt_pred <- predict(rf_15_wtd, newdata = form_test)
rf_wt_pred
cm_rf <- confusionMatrix(rf_wt_pred, as.factor(dat_test$clf))
cm_rf
### Try it with max depth 20 and no tune grid
# tunegrid <- expand.grid(.splitrule = "gini", .mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)

rf_20_wtd <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  verbose = TRUE,
  max.depth = 20,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_20_wtd)
saveRDS(rf_20_wtd, "data/20190305_rf_models/rf_20_wtd_small_sub_grp.rds")
             
rf_wt_pred <- predict(rf_20_wtd, newdata = form_test)
cm_rf <- confusionMatrix(rf_wt_pred, as.factor(dat_test$clf))
cm_rf
dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select[order(dtl_feat_select$F1),]

# write.csv(dtl_feat_select, "data/model_comparisons/rf_15_AA34_Precision_Recall_Accuracy_Functional_class_2012604_weighted.csv", row.names = T)


###########
# RF 20 unweighted
rf_20_unwtd <- caret::train(
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
  max.depth = 20,
  # class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_20_unwtd)

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
cm_rf$table

dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select[order(dtl_feat_select$F1),]
cm_df <- data.frame(melt(cm_rf$table), stringsAsFactors = F)
head(cm_df)

pdf("data/cm_rf_subst_grp_confusion_matrix.pdf")
ggplot(cm_df, aes(Reference, Prediction)) +
  geom_tile(aes(fill = value)) + 
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90))
dev.off()
