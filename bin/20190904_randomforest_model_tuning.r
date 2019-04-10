## Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/1742_aa34_signatures_20190104.fa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)

# source("src/convert_aln_15aap.r")
# aa <- convert_aln_15aap(rdaln) #5 physicochemical properties
# aadf <- data.frame(aa,stringsAsFactors = F)
# rownames(aadf)
# 
# # Write to CSV file
# write.csv(aadf, paste0("data/1742_seqs_510_feats_scaled_20190904.csv"), row.names=rownames(aadf), quote = F)

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



rf_grid <- expand.grid(mtry = c(2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
            splitrule = c("gini", "extratrees"),
            min.node.size = c(1, 3, 5))
            # max.depth = c(2, 3, 4, 5, 6, 7, 8, 9, 10))

# Tune grid for ranger
#### This is with an unlimited max depth
rf <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  tuneGrid = rf_grid,
  num.trees = 500,
  verbose = TRUE,
  importance = "permutation",
  preProcess = c("center", "scale"),
  classification = T)


saveRDS(rf, "data/20190904_model_comparisons/rf_xvalidated_tuned_no_max_depth_20190904.rds")


# Now with a max depth of 5 splits
rf5 <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  max.depth = 5,
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  tuneGrid = rf_grid,
  num.trees = 500,
  verbose = TRUE,
  importance = "permutation",
  preProcess = c("center", "scale"),
  classification = T)

saveRDS(rf5, "data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_5_20190904.rds")

# Now with a max depth of 7 splits
rf9 <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  max.depth = 9,
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  tuneGrid = rf_grid,
  num.trees = 500,
  verbose = TRUE,
  importance = "permutation",
  preProcess = c("center", "scale"),
  classification = T)

saveRDS(rf9, "data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_9_20190904.rds")


# Make a default model and test the effects of the number of trees on accuracy
rf_probs <- ranger(y_train ~., data = form_train, num.trees = 500, splitrule = "gini",
                 mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                 importance = "permutation", probability = TRUE)
saveRDS(rf_probs, "data/20190904_model_comparisons/rf_probability_ranger_20190904.rds")
rf_probs_pred <- predict(rf_probs, data = form_test)
rf_probs_pred$predictions
