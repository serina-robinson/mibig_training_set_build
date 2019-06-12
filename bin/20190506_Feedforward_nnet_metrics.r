# Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", 
               "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the data
rawdat <- read_csv("data/1553_training_sqs_with_loop_extracted.csv")
table(duplicated(rawdat[,2:ncol(rawdat)])) # check no duplicates
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[!grepl(paste0(c("HOLDOUT", "OTHER", "CAR", "amino.acid", "reject"), collapse = "|"), rawdat$clf),] # 658 observations
dat <- dat[,-grep("^X", colnames(dat))]
  

# SMOTE oversampling
# Stratified cross-validation

nnet_models <- list()
nnet_pred <- list()
dat_test <- list()
nnet_cm <- list()

for(i in 1:10) {
  
  set.seed(i)
  
  # Split into test and train
  dat_split <- initial_split(dat, prop = 3/4, strata = "clf")
  dat_train <- training(dat_split)
  dat_test  <- testing(dat_split)
  nrow(dat_train)/nrow(dat) # 75 %
  # table(word(dat$nms, sep = "_", 2))
  
  # Define our response
  x_train <- dat_train[,!colnames(dat_train) %in% c("nms", "clf")]
  x_test <- dat_test[,!colnames(dat_test) %in% c("nms", "clf")]
  y_train <- as.factor(dat_train$clf)
  y_test <- as.factor(dat_test$clf)

  print(table(y_test))
  
  # Complete dataset for training and testing
  form_train <- data.frame(cbind(x_train, y_train), stringsAsFactors = F, row.names = dat_train$nms)
  form_test <- data.frame(cbind(x_test, y_test), stringsAsFactors = F, row.names = dat_test$nms)
  
  classwts <- nrow(dat)/(length(unique(dat$clf)) * table(dat$clf))
  dtf <- data.frame(cbind(table(dat$clf), classwts))
  colnames(dtf) <- c("class_size", "class_weights")
  
  fitControl <- trainControl(method = "none", classProbs = FALSE)
  nnet_grid <- expand.grid(.decay = c(0.5),
                           .size = c(10))
  
  nnet_models[[i]] <- train(
    x = x_train,
    y = y_train,
    method = "nnet",
    MaxNWts = 100000,
    tuneGrid = nnet_grid,
    trControl = fitControl)
  
  nnet_models[[1]]
  nnet_pred[[i]] <- predict(nnet_models[[i]], x_test)
  nnet_cm[[i]] <- confusionMatrix(nnet_pred[[i]], as.factor(dat_test$clf))
  dat_test[[i]] <- dat_test
}

# Calculate the overall test statistics
length(nnet_cm)
xx <- list()
for (i in 1:length(nnet_cm)) {
  xx[[i]] <- nnet_cm[[i]]$overall[1]
}
mean(unlist(xx))
sd(unlist(xx))
