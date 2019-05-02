## Install packages
pacman::p_load("caret", "Biostrings", "RColorBrewer", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the data
rawdat <- read_csv("data/1742_seqs_510_feats_scaled_20190904.csv")
rawdat <- data.frame(cbind(rawdat$X1), scale(rawdat[,2:ncol(rawdat)]), stringsAsFactors = F)
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -3, sep = "_")
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid"), collapse = "|"), rawdat$nms),] # 658 observations

# Set seed 
set.seed(20190304)

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


# Make a default model and test the effects of the number of trees on accuracy
# rf_probs <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
#                    mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
#                    importance = "permutation", probability = TRUE)
# rf_grid <- expand.grid(mtry = c(2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
#                        splitrule = c("gini", "extratrees"),
#                        min.node.size = c(1, 3, 5))

rf <- caret::train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  num.trees = 500,
  verbose = TRUE,
  importance = "permutation")

rf$results

rf_probs <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                   mtry = 31, min.node.size = 1,
                   importance = "permutation", probability = TRUE)

saveRDS(rf, "data/rf_probability_small_subst_group_ntree1000_20191104.rds")


# Make a default model and test the effects of the number of trees on accuracy
rf_probs <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                   mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                   importance = "permutation", probability = TRUE)

# saveRDS(rf_probs, "data/rf_probability_ranger_ntree1000_20191104.rds")
rf_probs_pred <- predict(rf_probs, data = form_test)
rf_probs_pred$predictions
rowSums(rf_probs_pred$predictions)




colnames(rf_probs_pred$predictions)[apply(rf_probs_pred$predictions, 1, which.max)]

pred_df_test <- data.frame(cbind(as.character(y_test), colnames(rf_probs_pred$predictions)[apply(rf_probs_pred$predictions, 1, which.max)], apply(rf_probs_pred$predictions, 1, max)), stringsAsFactors = F)
head(pred_df_test)
# pred_df_train <- data.frame(cbind(as.character(y_train), colnames(rfprob$predictions)[apply(rfprob$predictions, 1, which.max)], apply(rfprob$predictions, 1, max)), stringsAsFactors = F)
colnames(pred_df_test) <- c("truth", "prediction", "probability")

pred_df_test$correct <- ifelse(pred_df_test$truth == pred_df_test$prediction, "Y", "N")
write.csv(pred_df_test, "false_testing_set_predictions.csv")
pred_df_test$probability <- as.numeric(pred_df_test$probability)
head(pred_df_test$probability)

length(pred_df_test$probability)
length(pred_df_test$probability[pred_df_test$probability > 0.5]) # Could try anything with a probability of prediction over 0.6 
344/425

summary(pred_df_test$probability[pred_df_test$correct == "Y"]) # Mean is 81%, median is 92
summary(pred_df_test$probability[pred_df_test$correct == "N"]) # Mean is 42% 

table(pred_df_test$correct)
