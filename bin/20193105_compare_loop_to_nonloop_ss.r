# Install packages
pacman::p_load("caret", "data.table", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20193105)

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
table(word(dat$nms, sep = "_", -1))
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
  max.depth = 20,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted)
rf_weighted$results
saveRDS(rf_weighted, "data/20193105_rf_weighted_with_loop_ss.rds")


### Try it with max depth 20 and tune grid determined
tunegrid <- expand.grid(.splitrule = "gini", .mtry = as.integer(sqrt(ncol(x_train))), .min.node.size = 1)

rf_unweighted <- caret::train(
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
  max.depth = 20,
  importance = "permutation")
getTrainPerf(rf_unweighted)
saveRDS(rf_unweighted, "data/20193105_rf_unweighted_with_loop_ss.rds")

##### Now try the exact same but modify to remove the extra amino acids from the loop sequence and compare performance
x_train <- x_train[,-grep("^X", colnames(x_train)),]
y_train <- y_train[,-grep("^X", colnames(y_train)),]

# Weighted
rf_weighted_no_loop <- caret::train(
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
  max.depth = 20,
  class.weights = classwts,
  importance = "permutation")
getTrainPerf(rf_weighted_no_loop)
saveRDS(rf_weighted_no_loop, "data/20193105_rf_weighted_no_loop_ss.rds")



# Unweighted
rf_unweighted_no_loop <- caret::train(
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
  max.depth = 20,
  importance = "permutation")
getTrainPerf(rf_unweighted_no_loop)
saveRDS(rf_unweighted_no_loop, "data/20193105_rf_unweighted_no_loop_ss.rds")


rf_weighted$results
rf_weighted_no_loop$results

rf_unweighted$results
rf_unweighted_no_loop$results

