# Install packages
pacman::p_load("caret", "reshape2", "RColorBrewer", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

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

# Read in the random forest models with different max depths
rfmax <- readRDS('data/20190904_model_comparisons/rf_xvalidated_tuned_no_max_depth_20190904.rds')
rf9 <- readRDS('data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_9_20190904.rds')
rf7 <- readRDS('data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_7_20190904.rds')
rf5 <- readRDS('data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_5_20190904.rds')
rf3 <- readRDS('data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_3_20190904.rds')
rf1 <- readRDS('data/20190904_model_comparisons/rf_xvalidated_tuned_max_depth_1_20190904.rds')
rfprob <- readRDS('data/20190904_model_comparisons/rf_probability_ranger_20190904.rds')
# All had different parameters, could try predicting with best parameter

# Rf prob
rfprob_pred <- predict(rfprob, data = form_test)

rowSums(rfprob$predictions)
dim(rfprob$predictions)
length(as.character(y_test))

colnames(rfprob$predictions)[apply(rfprob$predictions, 1, which.max)]

colnames(rfprob$predictions)[apply(rfprob$predictions, 1, which.max)]
dim(rfprob$predictions)

pred_df_train <- data.frame(cbind(as.character(y_train), colnames(rfprob$predictions)[apply(rfprob$predictions, 1, which.max)], apply(rfprob$predictions, 1, max)), stringsAsFactors = F)
colnames(pred_df_train) <- c("truth", "prediction", "probability")
pred_df_train$correct <- ifelse(pred_df_train$truth == pred_df_train$prediction, "Y", "N")
write.csv(pred_df_train, "false_training_set_predictions.csv")
pred_df_train$probability <- as.numeric(pred_df_train$probability)
head(pred_df_train$probability)

length(pred_df_train$probability[pred_df_train$probability > 0.6]) # Could try anything with a probability of prediction over 0.6 

summary(pred_df_train$probability[pred_df_train$correct == "Y"]) # Mean is 81%, median is 92
summary(pred_df_train$probability[pred_df_train$correct == "N"]) # Mean is 42% 

table(pred_df_train$correct)
# cm_rfmax <- confusionMatrix(rfprob_pred$predictions, as.factor(dat_test$clf))

# Rf max
rfmax_pred <- predict(rfmax, newdata = form_test)
cm_rfmax <- confusionMatrix(rfmax_pred, as.factor(dat_test$clf))

# Rf 9
rf9_pred <- predict(rf9, newdata = form_test)
cm_rf9 <- confusionMatrix(rf9_pred, as.factor(dat_test$clf))

# Rf 7
rf7_pred <- predict(rf7, newdata = form_test)
cm_rf7 <- confusionMatrix(rf7_pred, as.factor(dat_test$clf))

# Rf 5
rf5_pred <- predict(rf5, newdata = form_test)
cm_rf5 <- confusionMatrix(rf5_pred, as.factor(dat_test$clf))

# Rf 3
rf3_pred <- predict(rf3, newdata = form_test)
cm_rf3 <- confusionMatrix(rf3_pred, as.factor(dat_test$clf))

# Rf 1
rf1_pred <- predict(rf1, newdata = form_test)
cm_rf1 <- confusionMatrix(rf1_pred, as.factor(dat_test$clf))

mylist <- list(cm_rfmax$byClass, cm_rf9$byClass, cm_rf7$byClass, cm_rf5$byClass, cm_rf3$byClass, cm_rf1$byClass)

F1 <- do.call("cbind", (lapply(1:length(mylist), function(x) { mylist[[x]][,"F1"]})))
bal_acc <- do.call("cbind", (lapply(1:length(mylist), function(x) { mylist[[x]][,"Balanced Accuracy"]})))
colnames(bal_acc) <- c("No_limit", "Limit_depth_9", "Limit_depth_7", "Limit_depth_5", "Limit_depth_3", "Limit_depth_1")
rownames(bal_acc) <- gsub("Class: ", "", rownames(bal_acc))

bal_long <- reshape2::melt(bal_acc)
head(bal_long)
dim(bal_long)
head(bal_long)

pal1 <- colorRampPalette(colors=brewer.pal(8, "Set1"))(8)
pal1
pal1[pal1 == "#984EA3"] <- "firebrick1"
pal1[2] <- "#984EA3"
pal1[pal1 == "#E41A1C"] <- "#377EB8"
pal1[pal1 == "#A65628"] <- "chartreuse3"
# pa1l1[pal1 == ""] 

pal1[pal1 == "#FFFF33"] <- "goldenrod4"
pal2 <- c(pal1, "navy", #"black", 
          "gray68", "plum1", "blue1",
          "deepskyblue", "gold", "darkorchid1", 
          "deeppink2", "lightslateblue",
          "lightblue2", "darkseagreen1")

palette(pal2)
pdf("output/impact_of_tree_depth.pdf")
ggplot(bal_long, aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_point() +
  geom_path() +
  ylab("Within group classification accuracy") +
  xlab("Maximum tree depth limit") +
  scale_color_manual(values = pal2) +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 90),
        legend.text = element_text(size = 8),
        legend.title = element_blank())
dev.off()


# Look at other parameters
rfmax$control
rfmax$results$AccuracySD
rfmax_res <- rfmax$results %>%
  dplyr::select(mtry, splitrule, min.node.size, Accuracy, AccuracySD)

findf <- rfmax_res %>%
  dplyr::filter(splitrule == "gini") %>%
  dplyr::filter(min.node.size == 1)
findf


pdf("output/impact_of_mtry_on_accuracy.pdf")
ggplot(findf, aes(x = mtry, y = Accuracy, color = mtry)) +
  geom_point() +
  xlab("Number of variables available for splitting at each tree node") +
  ylab("Within class accuracy (test set)") +
  geom_errorbar(aes(ymax = Accuracy + AccuracySD, ymin = Accuracy - AccuracySD)) +
  ylim(0.1, 1)
dev.off()
