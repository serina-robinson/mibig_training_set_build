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
head(rawdat)
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)

# Try downsampling the ARYL and NRPS sequecnes to about 80 
# 4-coumarate sequences in the aryl class
# sample(which(rawdat$clf == "NRPS"), 80)
# sample(which(rawdat$clf == "ARYL"), 80)

# rawdat_tr <- rawdat %>%
#   dplyr::filter(stringr::str_detect(pattern = "ARYL|NRPS", clf, negate = T)) %>%
#   bind_rows(., rawdat[sample(which(rawdat$clf == "NRPS"), 80),]) %>%
#   bind_rows(., rawdat[sample(which(rawdat$clf == "ARYL"), 80),])
# dim(rawdat_tr)

# Remove the holdout test predictions
# dat <- rawdat_tr[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat_tr$nms),] # 658 observations
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat$nms),]
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
saveRDS(rf, "data/20190305_rf_models/20190305_rf_no_wts_no_max_depth_funct_class.rds")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
cm_rf

sink("data/model_comparisons/cm_rf_functional_class_20190503.txt")
cm_rf
sink()

dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_Precision_Recall_Accuracy_Functional_class_20190503.csv", row.names = T)

# Try with weights
classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
classwts
dtf <- data.frame(cbind(table(dat$clf), classwts))
colnames(dtf) <- c("class_size", "class_weights")
# write.csv(dtf, "output/small_sub_grp_class_sizes_and_weights.csv", quote = F, row.names = T)

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
  class.weights = classwts,
  # max.depth = 20,
  importance = "permutation")
getTrainPerf(rf)


# saveRDS(rf, "data/20190305_rf_models/rf_wtd_no_max_depth_funct_class.rds")
rf_pred2 <- predict(rf, newdata = form_test)
cm_rf2 <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
cm_rf2$overall
cm_rf$overall
rf$bestTune

tunegrid <- expand.grid(.mtry =  34, .splitrule = 'gini', .min.node.size = 1)
# Now try with max depth 20
rf_20_wtd <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 3,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 1000,
  tuneGrid = tunegrid,
   #as.integer(sqrt(ncol(x_train))), .min.node.size = 1),
  verbose = TRUE,
  max.depth = 20,
  class.weights = classwts,
  importance = "permutation")

# Max depth 20 weighted
saveRDS(rf_20_wtd, "data/20190305_rf_models/rf_wtd_depth_20.rds")


rf_wt_pred <- predict(rf_20_wtd, newdata = form_test)
cm_rf <- confusionMatrix(rf_wt_pred, as.factor(dat_test$clf))
cm_rf

sink("data/model_comparisons/cm_rf_wt_max_depth_20_pred_20190305.txt")
cm_rf
sink()

dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "F1", "Balanced Accuracy")], 2))
dtl_feat_select
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_Precision_Recall_Accuracy_Functional_class_20190305_weighted.csv", row.names = T)


