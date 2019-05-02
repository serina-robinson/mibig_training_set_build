## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20191304)

# Read in the data
rawdat <- read_csv("data/797_seqs_510_feats_for_supervised_scaled_20191104.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")
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

rf_fc <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation", probability = TRUE) 


saveRDS(rf_fc, "data/rf_functional_class_1000trees_probability.rds")

rf_fc_no_prob <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                        mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                        importance = "permutation")
rf_fc_pred <- predict(rf_fc_no_prob, data = form_test)

cm_rf_fc <- confusionMatrix(rf_fc_pred$predictions, as.factor(dat_test$clf))
cm_rf_fc$byClass

######### Now by substrate group #########

# Read in the data
rawdat <- read_csv("data/1742_seqs_510_feats_scaled_20190904.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -3, sep = "_")
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

rf_fc <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation", probability = TRUE) 


saveRDS(rf_fc, "data/rf_substrate_group_1000trees_probability.rds")

rf_fc_no_prob <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                        mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                        importance = "permutation")
rf_fc_pred <- predict(rf_fc_no_prob, data = form_test)

cm_rf_fc <- confusionMatrix(rf_fc_pred$predictions, as.factor(dat_test$clf))
cm_rf_fc$byClass
