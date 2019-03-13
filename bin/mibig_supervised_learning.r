## Install packages
pacman::p_load("caret", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")
# devtools::install_github('mlampros/FeatureSelection')
# devtools::install_version("xgboost", version = "0.4-4", repos = "http://cran.us.r-project.org")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

##### Reading in the data #####
rawdat <- read_csv("data/710_seqs_1592_feats_for_supervised.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")

# Consolidate some classes that are too small
clndat <- rawdat %>%
  mutate(nms_cln = case_when(clf == "NOTE" ~ paste0(nms, "_NRPS"),
                             clf == "COUM" ~ paste0(nms, "_ARYL"),
                             clf == "BIARYL" ~ paste0(nms, "_ARYL"),
                             clf == "MMCS" ~ paste0(nms, "_SACS"),
                             TRUE ~ nms))

rawdat$clf <- word(clndat$nms_cln, -1, sep = "_") 
table(rawdat$clf)

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


#### Random Forest without feature selection ####

# rf <- train(
#   x = x_train,
#   y = y_train,
#   method = "ranger",
#   trControl = trainControl(method = "cv", number = 10,
#                            verboseIter = T, classProbs = T,
#                            savePredictions = "final"),
#   # tuneGrid = tgrid,
#   num.trees = 500,
#   # respect.unordered.factors = FALSE,
#   verbose = TRUE,
#   importance = "permutation")


# Aggregating results
# Selecting tuning parameters
# Fitting mtry = 56, splitrule = gini, min.node.size = 1 on full training set

# Ranger 
rf_no_lasso <- ranger(classification = T, dependent.variable.name = 'y_train',
                data = form_train, num.trees = 500, mtry = 56,
                splitrule = "gini", min.node.size = 1,
                importance = "permutation")

rf_pred <- predict(rf_no_lasso, data = form_test)
rf_pred$predictions

cm_rf_no_lasso <- confusionMatrix(rf_pred$predictions, form_test$y_test)
sink("data/cm_rf_no_lasso.txt")
print(cm_rf_no_lasso)
sink()

# Variable importance plot
rf_no_lasso$variable.importance

rf_imp <- sort(rf_no_lasso$variable.importance, decreasing = T)
rf_imp
# names(rf_imp) <- factor(names(rf_imp), levels = names(rf_imp)[sort(rf_imp, decreasing = T)])
# rf_imp[1:12] <- factor(rf_imp[1:12], levels = rf_imp[1:12])
# rf_imp[1:12]
vimp <- data.frame(cbind(sort(rf_no_lasso$variable.importance, decreasing = T), 
                         names(sort(rf_no_lasso$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
         variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp
pdf("data/rf_no_lasso_var_imp.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
ggsave("data/rf_no_lasso_var_imp.jpeg")
dev.off()



# Feature selection using LASSO
library(doMC)
registerDoMC(2)
x_cent_scale <- as.matrix(x_train)

# LASSO for feature selection
cv1 <- cv.glmnet(x_cent_scale, y_train, family = "multinomial", nfold = 10, parallel = TRUE, alpha = 1)
keep <- coef(cv1, cv1$lambda.1se)
pos_coefs <- list()

for(i in 1:length(keep)) {
  tmp <- as.matrix(keep[[i]])
  pos_coefs[[i]] <- rownames(tmp)[tmp != 0]
}

feats_with_intercept <- unique(unlist(pos_coefs))
final_feats <- feats_with_intercept[feats_with_intercept != "(Intercept)"]
final_feats # 183 feats instead of 1592

rownames(x_cent_scale) <- rownames(df_train)
feat_select_train <- x_cent_scale[,colnames(x_cent_scale) %in% final_feats]

#### Random Forest with feature selection ####
# Random Forest using the ranger package for prediction
# Fitting mtry = 2, splitrule = gini, min.node.size = 1 on full training set
# rf <- train(
#   x = feat_select_train, 
#   y = y_train,
#   method = "ranger",
#   trControl = trainControl(method = "cv", number = 10, 
#                            verboseIter = T, classProbs = T,
#                            savePredictions = "final"),
#   # tuneGrid = tgrid,
#   num.trees = 500,
#   # respect.unordered.factors = FALSE,
#   verbose = TRUE,
#   importance = "permutation")
# 
# 
# Variable importance plot
rf_imp <- varImp(rf, scale = FALSE,
                 surrogates = FALSE,
                 competes = FALSE)

pdf("data/rf_var_imp.pdf", width = 3, height = 3)
ggplot(rf_imp, top = 12) + xlab("") + theme_bw()
dev.off()
# 
# # ROC curve
# approx_roc_curve <- function(x, label) {
#   x %>%
#     pluck("pred") %>%
#     roc_curve(obs, y_train) %>%
#     mutate(model = label)
# }
# 
# # ROC curve is not looking good...
# approx_roc_curve(rf, "Random Forest") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)
# 
# # Confusion matrix
# getTrainPerf(rf)
# 

# Select the same features in the test set
dim(dat_test)
feat_select_test <- dat_test[,colnames(dat_test) %in% colnames(feat_select_train)]
feat_select_test_final <- data.frame(cbind(feat_select_test, dat_test$clf))
dim(feat_select_test_final)

# Prediction
class(feat_select_train)
feat_select_train_full <- data.frame(feat_select_train) %>%
  tibble::add_column(., y_train)
feat_select_train_full$y_train
rf_lasso_select <- ranger(classification = T, dependent.variable.name = 'y_train',
                data = feat_select_train_full, num.trees = 500, mtry = 2,
                splitrule = "gini", min.node.size = 1,
                importance = "permutation")
rf_pred <- predict(rf_lasso_select, data = feat_select_test_final)

cm_rf <- confusionMatrix(rf_pred$predictions, as.factor(dat_test$clf))
cm_rf
sink("data/cm_rf_feat_select.txt")
print(cm_rf)
sink()


vimp <- data.frame(cbind(sort(rf_lasso_select$variable.importance, decreasing = T), 
                         names(sort(rf_lasso_select$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
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


dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "output/Feat_Select_Precision_Recall_Accuracy.csv", row.names = T)

dtl_no_lasso <- data.frame(round(cm_rf_no_lasso$byClass[,colnames(cm_rf_no_lasso$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_no_lasso, "output/No_Lasso_Precision_Recall_Accuracy.csv", row.names = T)
dtl_no_lasso

## Make a heatmap of confusion matrix results

cm_list <- list(cm_rf$table,
                cm_rf_no_lasso$table)
names(cm_list) <- c("rf_feat_selection_lasso", "rf_all_feats")

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
plot_grid(pllist[[1]],
          pllist[[2]])
ggsave("conf_matrices.jpeg", height=4, width=10, units='in')
# dev.off()
# ggsave("conf_matrices.jpeg", height=6, width=10, units='in')



rf_pred <- predict(rf, data = feat_select_test_final)dim(feat_select_train)  