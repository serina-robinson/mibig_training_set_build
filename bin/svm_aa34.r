## Install packages
pacman::p_load("caret", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the data
rawdat <- read_csv("data/713_seqs_136_feats_for_supervised.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER"), collapse = "|"), rawdat$nms),] # 658 observations

# Set seed 
set.seed(20190503)

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

#### SVM with radial basis kernel ####
# Fitting sigma = 0.00396, C = 0.25 on full training set
svm_radial <- train(
  x = df_train,
  y = y_train,
  method = "svmRadial",
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  preProcess = c("center", "scale"),
  verbose = TRUE)


# Variable importance plot
svm_radial_imp <- varImp(svm_radial, scale = FALSE,
                  surrogates = FALSE,
                  competes = FALSE)


pdf("data/svm_var_imp.pdf", width = 10, height = 10)
ggplot(svm_imp, top = 12) + xlab("") + theme_bw()
dev.off()


# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}

# ROC curve is not looking good...
approx_roc_curve(svm_radial, "SVM Radial Kernel") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(svm_radial)

# Select the same features in the test set

# Try prediction
svm_ml <- svm(y_train ~., data = form_train, kernel = "radial",
              sigma = 0.00396, C = 0.25)
summary(svm_ml)

svm_pred <- predict(svm_ml, newdata = form_test)
length(svm_pred)

cm_svm_radial <- confusionMatrix(svm_pred, as.factor(dat_test$clf))
cm_svm_radial

sink("data/cm_svm_feat_select.txt")
print(cm_svm_radial)
sink()

dtl_feat_select <- data.frame(round(cm_svm_radial$byClass[,colnames(cm_svm_radial$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "output/AA34_SVM_Radial_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)


## Make a heatmap of confusion matrix results

cm_list <- list(cm_svm_radial$table)
names(cm_list) <- c("svm_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
  #axis.title.y = element_blank())
}

# jpeg("output/data_confusion_matrices.jpeg", width = 750, height = 500, res = 300)
pllist[[1]]
ggsave("aa34_group_svm_conf_matrices.jpeg", height=4, width=5, units='in')
