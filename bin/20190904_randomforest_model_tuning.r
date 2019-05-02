## Install packages
pacman::p_load("caret", "Biostrings", "RColorBrewer", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

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



# Now test the effect of number of trees in the case where 
for(i in 1:10) {
  rf1500 <- ranger(y_train ~., data = form_train, num.trees = 1500, splitrule = "gini",
          mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
          importance = "permutation") 
  saveRDS(rf1500, paste0("data/20191104_ntree_comparisons/rf_1500_trees_fold_", i, ".rds"))
  rf1000 <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation")
  saveRDS(rf1000, paste0("data/20191104_ntree_comparisons/rf_1000_trees_fold_", i, ".rds"))
  rf750 <- ranger(y_train ~., data = form_train, num.trees = 750, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation")
  saveRDS(rf750, paste0("data/20191104_ntree_comparisons/rf_750_trees.rds_fold_", i, ".rds"))
  rf500 <- ranger(y_train ~., data = form_train, num.trees = 500, splitrule = "gini",
                   mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                   importance = "permutation")
  saveRDS(rf500, paste0("data/20191104_ntree_comparisons/rf_500_trees.rds_fold_", i, ".rds"))
  
  rf250 <- ranger(y_train ~., data = form_train, num.trees = 100, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation")
  saveRDS(rf250, paste0("data/20191104_ntree_comparisons/rf_250_trees.rds_fold_", i, ".rds"))
  
  rf100 <- ranger(y_train ~., data = form_train, num.trees = 100, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation")
  saveRDS(rf100, paste0("data/20191104_ntree_comparisons/rf_100_trees_fold_", i, ".rds"))
  
  rf50 <- ranger(y_train ~., data = form_train, num.trees = 50, splitrule = "gini",
                mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                importance = "permutation")
  saveRDS(rf50, paste0("data/20191104_ntree_comparisons/rf_50_trees_fold_", i, ".rds"))
  
  rf10 <- ranger(y_train ~., data = form_train, num.trees = 10, splitrule = "gini",
               mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
               importance = "permutation")
  saveRDS(rf10, paste0("data/20191104_ntree_comparisons/rf_10_trees_fold_", i, ".rds"))
}

# saveRDS(rf100, "data/20190904_model_comparisons/rf_100_trees.rds")

# Rf 
rffils <- list.files("data/20191104_ntree_comparisons/", full.names = T)
length(rffils)

rftrees <- list()
rffils
mods <- lapply(1:length(rffils), function(x) {readRDS(rffils[x])})
mods[[80]]$prediction.error

# mods <- rftrees
rftrees <- lapply(1:length(mods), function(x) { predict(mods[[x]], data = form_test) })
mylist <- lapply(1:length(mods), function(x) confusionMatrix(rftrees[[x]]$predictions, as.factor(dat_test$clf)))

# mylist <- list(cm_rfmax$byClass, cm_rf9$byClass, cm_rf7$byClass, cm_rf5$byClass, cm_rf3$byClass, cm_rf1$byClass)

F1 <- do.call("cbind", (lapply(1:length(mylist), function(x) { mylist[[x]]$byClass[,"F1"]})))
# bal_acc <- do.call("cbind", (lapply(1:length(mylist), function(x) { mylist[[x]]$byClass[,"Balanced Accuracy"]})))
dim(bal_acc) # 15 rows by 80 columns

bal_acc <- do.call("cbind", (lapply(1:length(mods), function(x) { mods[[x]]$prediction.error })))
bal_acc
# colnames(bal_acc) <- rep(c(1500, 1000, 750, 500, 250, 100, 50, 10), each = 10)
colnames(bal_acc) <- rep(c(10, 100, 1000, 1500, 250, 50, 500, 750), each = 10)
colnames(bal_acc)
# rownames(bal_acc) <- gsub("Class: ", "", rownames(bal_acc))

bal_long <- reshape2::melt(bal_acc)
head(bal_long)
bal_long$Var1
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

head(bal_long)
bal_mean <- bal_long %>%
  # group_by(Var1, Var2) %>%
  group_by(Var2) %>%
  summarize(mean_acc = mean(value, na.rm = TRUE))
bal_sd <- bal_long %>%
  group_by(Var2) %>%
  summarize(mean_sd = sd(value, na.rm = TRUE))

# limits <- c(ymin = bal_mean$mean_acc - bal_sd$mean_sd, ymax = bal_mean$mean_acc + bal_sd$mean_sd )
# limits

pdf("output/impact_of_num_trees_pred_error.pdf", width = 5, height = 5)
ggplot(bal_mean, aes(x = Var2, y = mean_acc)) + #color = Var1, group = Var1)) +
  geom_point() +
  geom_path() +
  geom_errorbar(aes(ymin = bal_mean$mean_acc - bal_sd$mean_sd, ymax = bal_mean$mean_acc + bal_sd$mean_sd )) +
  ylab("OOB Prediction Error") +
  xlab("Number of trees") +
  scale_color_manual(values = pal2) +
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 90),
        legend.text = element_text(size = 8),
        legend.title = element_blank()) +
  ylim(0.15, 0.25)
dev.off()





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
rf_probs <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                 mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                 importance = "permutation", probability = TRUE)
saveRDS(rf_probs, "data/20190904_model_comparisons/rf_probability_ranger_20190904.rds")
rf_probs_pred <- predict(rf_probs, data = form_test)
rf_probs_pred$predictions
