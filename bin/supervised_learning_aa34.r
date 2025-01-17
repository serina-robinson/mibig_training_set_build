## Install packages
pacman::p_load("caret", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the data
rawdat <- read_csv("data/797_seqs_136_feats_for_supervised_20190328.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -3, sep = "_")

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("HOLDOUT", "OTHER", "CAR"), collapse = "|"), rawdat$nms),] # 658 observations

# Set seed 
set.seed(20192603)

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

#### Random Forest with feature selection ####
# Random Forest using the ranger package for prediction
# Fitting mtry = 2, splitrule = gini, min.node.size = 1 on full training set
rf <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation")

# 
# Variable importance plot
rf_imp <- varImp(rf, scale = FALSE,
                 surrogates = FALSE,
                 competes = FALSE)

pdf("data/rf_var_imp.pdf", width = 3, height = 3)
ggplot(rf_imp, top = 12) + xlab("") + theme_bw()
dev.off()
# 
# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}

# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# Select the same features in the test set

# Try prediction


rf_ml <- ranger(y_train ~., data = form_train, num.trees = 500, splitrule = as.character(rf$bestTune$splitrule), 
                mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")
# saveRDS(rf_ml, "output/rf_functional_class_20192803.rds")
rf_ml <- readRDS("output/rf_small_subst_grp_20192803.rds")

alldat <- data.frame(bind_rows(form_train, form_test))
actual_clf <- c(as.character(y_train), as.character(y_test))
nms <- c(dat_train$nms, dat_test$nms)
nms

rf_pred <- predict(rf_ml, data = alldat)
length(rf_pred$predictions)


cm_rf <- confusionMatrix(rf_pred$predictions, as.factor(actual_clf))
cm_rf

sink("data/cm_rf_feat_select.txt")
print(cm_rf)
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T), 
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp

pdf("data/rf_feat_select_var_imp.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
ggsave("data/rf_feat_select_var_imp.jpeg")
dev.off()


# dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
# write.csv(dtl_feat_select, "output/AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

dtf <- data.frame(cbind(nms, as.character(rf_pred$predictions), as.character(actual_clf)), stringsAsFactors = F)
dtf
colnames(dtf) <- c("sqnams_tr", "pred", "actual")
dtf_wrong <- dtf[dtf$pred != dtf$actual,]
dtf_wrong
ad <- read_csv("data/combined_adenylate_forming_training_set_20192703.csv")
dtf_wrong

ad_joined <- inner_join(dtf_wrong, ad, by = "sqnams_tr")
dim(ad_joined)

#ad_wrong <- data.frame(cbind(dtf_wrong[1:37,], ad[ad$sqnams_tr %in% dtf_wrong$nms,]))
write_csv(ad_joined, "data/classification_group_predictions_20192803.csv")

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
  #axis.title.y = element_blank())
}


# jpeg("output/data_confusion_matrices.jpeg", width = 750, height = 500, res = 300)
pllist[[1]]
ggsave("aa34_group_conf_matrices.jpeg", height=10, width=10, units='in')
