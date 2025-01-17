###Naive Bayes method
nb_grid <- expand.grid(usekernel = TRUE, fL = 0, adjust = 1)

nb <- train(
  x = df_train, 
  y = y_train,
  method = "nb",
  tuneGrid = nb_grid,
  trControl = trainControl(method = "repeatedcv", number = 10, 
                           repeats = 5,
                           verboseIter = T, classProbs = F,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nb_ml)
nb_ml$results

saveRDS(nb, "data/model_comparisons/nb_NRPS_broad_monomers_repeatedCV_20190306.rds")
nb_ml <- readRDS("data/model_comparisons/nb_NRPS_broad_monomers_repeatedCV_20190306.rds")
nb_pred <- predict(nb_ml, newdata = form_test)

cm_nb <- confusionMatrix(nb_pred, as.factor(dat_test$clf))
cm_nb

dtl_feat_select <- data.frame(round(cm_nb$byClass[,colnames(cm_nb$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/nb_AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nb$table)
names(cm_list) <- c("nb_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/model_comparisons/aa34_nb_conf_matrices_subst_grp_repeatedCV_20190306.jpeg", height=7, width=7, units='in')


####### Random Forest with feature selection ########
# Random Forest using the ranger package for prediction
# Fitting mtry = 2, splitrule = gini, min.node.size = 1 on full training set
rf <- train(
  x = x_train,
  y = y_train,
  method = "ranger",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation")
# preProcess = c("center", "scale"))

saveRDS(rf, "data/model_comparisons/rf_xvalidated_small_subst_grp_repeatedCV_20190306.rds")

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

# ROC curve
approx_roc_curve <- function(x, label) {
  x %>%
    pluck("pred") %>%
    roc_curve(obs, y_train) %>%
    mutate(model = label)
}
# saveRDS(rf, "data/model_comparisons/rf_xvalidated_small_subst_grp_20190306.rds")

# ROC curve is not looking good...
approx_roc_curve(rf, "Random Forest") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

# Confusion matrix
getTrainPerf(rf)

# Try prediction
rf_ml <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = as.character(rf$bestTune$splitrule), 
                mtry = rf$bestTune$mtry, min.node.size = rf$bestTune$min.node.size,
                importance = "permutation")
saveRDS(rf_ml, "data/model_comparisons/rf_NRPS_broad_monomers_repeatedCV_20190306.rds")
rf_ml <- readRDS("data/model_comparisons/rf_NRPS_broad_monomers_repeatedCV_20190306.rds")

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))

sink("data/model_comparisons/cm_rf_15aa_functional_class_repeatedCV_20190306.txt")
cm_rf
sink()

vimp <- data.frame(cbind(sort(rf_ml$variable.importance, decreasing = T), 
                         names(sort(rf_ml$variable.importance, decreasing = T))), stringsAsFactors = F) %>%
  dplyr::rename(importance = X1,
                variable = X2) %>%
  mutate(importance = as.numeric(importance)) %>%
  dplyr::slice(1:25)
vimp

pdf("data/model_comparisons/rf_feat_select_var_imp_NRPS_broad_monomers_repeatedCV_20190306.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()


dtl_feat_select <- data.frame(round(cm_rf$byClass[,colnames(cm_rf$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/rf_AA34_NRPS_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

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
}

pllist[[1]]
ggsave("data/model_comparisons/aa34_rf_NRPS_conf_matrices_subst_grp_20190306.jpeg", height=7, width=7, units='in')


###==================================================================###
###Neural networks using nnet
# Currently throws an error
nnet_grid <- expand.grid(.decay = c(0.5, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7),
                         .size = c(3, 5, 10, 20))

nnet_mod <- train(
  x = df_train,
  y = y_train,
  method = "nnet",
  MaxNWts = 10000,
  # size = 1,
  tuneGrid = nnet_grid,
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T,
                           savePredictions = "final"))

# Confusion matrix
getTrainPerf(nnet_mod)
# confusionMatrix(nnet_mod)

saveRDS(nnet_mod, "data/model_comparisons/nnet_NRPS_broad_monomers_repeatedCV_20190306.rds")
nnet_ml <- readRDS("data/model_comparisons/nnet_NRPS_broad_monomers_repeatedCV_repeatedCV_20190306.rds")
nnet_ml$bestTune
nnet_ml$call
# size decay
#  10 0.5

# approx_roc_curve(nnet_ml, "Neural network") %>%
#   ggplot(aes(x = 1 - specificity, y = sensitivity)) +
#   geom_path()  +
#   geom_abline(col = "red", alpha = .5)

nnet_pred <- predict(nnet_ml, newdata = form_test)
cm_nnet <- confusionMatrix(nnet_pred, as.factor(dat_test$clf))
cm_nnet

sink("data/model_comparisons/cm_nnet_NRPS_15aa_functional_class_20190306.txt")
cm_nnet
sink()

dtl_feat_select <- data.frame(round(cm_nnet$byClass[,colnames(cm_nnet$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/nnet_NRPS_AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_nnet$table)
names(cm_list) <- c("nnet_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/model_comparisons/aa34_nnet_conf_matrices_subst_grp_20190306.jpeg", height=7, width=7, units='in')


###==================================================================###
###SVM Radial basis kernel
# Fitting sigma = 0.000646, C = 1 on full training set
grid <- expand.grid(.C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
                    .sigma = c(0.1, 0.25,0.5,0.75,1))

svm_radial <- train(
  x = df_train,
  y = y_train,
  method = 'svmRadial',
  tuneGrid = grid,
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = F,
                           savePredictions = "final"),
  verbose = TRUE)

# Confusion matrix
getTrainPerf(svm_radial)
svm_radial$bestTune
saveRDS(svm_radial, "data/model_comparisons/svm_radial_xvalidated_small_subst_grp_20190306.rds")

approx_roc_curve(svm_radial, "SVM Radial") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

pdf("data/model_comparisons/svm_radial_feat_select_var_imp_NRPS_broad_monomers_repeatedCV_repeatedCV_20190306.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()

cm_svm_radial_pred <- predict(svm_radial, newdata = form_test)
cm_svm_radial <- confusionMatrix(svm_radial_pred, as.factor(dat_test$clf))

sink("data/model_comparisons/cm_nnet_15aa_functional_class_20190306.txt")
cm_nnet
sink()
dtl_feat_select <- data.frame(round(cm_svm_radial$byClass[,colnames(cm_svm_radial$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/svm_radial_AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_svm_radial$table)
names(cm_list) <- c("svm_radial_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/aa34_svm_radial_conf_NRPS_matrices_subst_grp_20190306.jpeg", height=7, width=7, units='in')

###==================================================================###
###SVM Linear kernel
# Fitting sigma = 0.000646, C = 1 on full training set


grid <- expand.grid(C = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
                    gamma=c(0.1, 0.25,0.5,0.75,1))
svm_linear <- train(
  x = df_train,
  y = y_train,
  tuneGrid= grid,
  tuneLength = 10,
  method = "svmLinear",
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = F,
                           savePredictions = "final")
  verbose = TRUE)

# Confusion matrix
getTrainPerf(svm_linear)
saveRDS(svm_linear, "data/model_comparisons/svm_linear_xvalidated_small_subst_grp_20190306.rds")

approx_roc_curve(svm_linear, "SVM Linear") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

pdf("data/model_comparisons/svm_radial_feat_select_var_imp_NRPS_broad_monomers_repeatedCV_20190306.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()

cm_svm_radial_pred <- predict(svm_radial, newdata = form_test)
cm_svm_radial <- confusionMatrix(svm_radial_pred, as.factor(dat_test$clf))

sink("data/model_comparisons/cm_nnet_15aa_functional_class_20190306.txt")
cm_nnet
sink()
dtl_feat_select <- data.frame(round(cm_svm_radial$byClass[,colnames(cm_svm_radial$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/svm_radial_AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_svm_radial$table)
names(cm_list) <- c("svm_radial_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/aa34_svm_radial_conf_matrices_subst_grp_20190306.jpeg", height=7, width=7, units='in')

##  Compare all confusion matrices
cm_list <- list(cm_rf$table,
                cm_nb$table,
                cm_nnet$table,
                cm_svm_radial$table)
names(cm_list) <- c("rf", "nb", "nnet", "svm_radial")

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

pllist[[1]]
# pdf("output/model_confusion_matrices.pdf", width = 13, height = 10)
jpeg("output/model_confusion_matrices.jpeg", width = 750, height = 500, res = 300)
plot_grid(pllist[[1]],
          pllist[[2]],
          pllist[[3]],
          pllist[[4]])
ggsave("conf_matrices.jpeg", height=6, width=10, units='in')

##  Compare all confusion matrices
cm_list <- list(cm_rf$table,
                cm_nb$table,
                cm_nnet$table,
                cm_svm_radial$table)
names(cm_list) <- c("rf", "nb", "nnet", "svm_radial")

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

pllist[[1]]
# pdf("output/model_confusion_matrices.pdf", width = 13, height = 10)
jpeg("output/model_confusion_matrices.jpeg", width = 750, height = 500, res = 300)
plot_grid(pllist[[1]],
          pllist[[2]],
          pllist[[3]],
          pllist[[4]])
ggsave("conf_matrices.jpeg", height=6, width=10, units='in')


###==================================================================###
###SVM Linear basis kernel
# Fitting sigma = 0.000646, C = 1 on full training set
grid <- expand.grid(.C = c(0.05, 0.1, 0.25, 0.5, 0.75, 1))
#sigma = c(0.1, 0.25,0.5,0.75,1))
grid

svm_linear <- train(
  x = df_train,
  y = y_train,
  method = "svmLinear",
  tuneGrid = grid,
  trControl = trainControl(method = "repeatedcv", number = 10,
                           repeats = 5,
                           verboseIter = T, classProbs = F,
                           savePredictions = "final")
  verbose = TRUE)

# Confusion matrix
getTrainPerf(svm_linear)
saveRDS(svm_linear, "data/model_comparisons/svm_linear_xvalidated_small_subst_grp_20190306.rds")

approx_roc_curve(svm_linear, "SVM Linear") %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path()  +
  geom_abline(col = "red", alpha = .5)

pdf("data/model_comparisons/svm_linear_feat_select_var_imp_NRPS_broad_monomers_repeatedCV_20190306.pdf", width = 6, height = 6)
ggplot(data = vimp, 
       aes(x=reorder(variable,importance), y=importance, fill=importance))+ 
  geom_bar(stat="identity", position="dodge")+ coord_flip()+
  ylab("Variable Importance")+
  xlab("")+
  guides(fill=F)+
  scale_fill_gradient(low="red", high="blue")
dev.off()

cm_svm_linear_pred <- predict(svm_linear, newdata = form_test)
cm_svm_radial <- confusionMatrix(svm_radial_pred, as.factor(dat_test$clf))

sink("data/model_comparisons/cm_nnet_15aa_monomers_20190306.txt")
cm_nnet
sink()
dtl_feat_select <- data.frame(round(cm_svm_radial$byClass[,colnames(cm_svm_radial$byClass) %in% c("Precision", "Recall", "Balanced Accuracy")], 2))
write.csv(dtl_feat_select, "data/model_comparisons/svm_radial_AA34_Precision_Recall_Accuracy_Substrate_group.csv", row.names = T)

## Make a heatmap of confusion matrix results
cm_list <- list(cm_svm_linear$table)
names(cm_list) <- c("svm_linear_aa_34")

pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[i]]), aes(Prediction, Reference)) +
    geom_tile(aes(fill = Freq)) + 
    geom_text(aes(label = round(Freq, 1))) +
    scale_fill_gradient(low = "white", high = "red") +
    ggtitle(names(cm_list)[i]) +
    theme(axis.text.x = element_text(angle = 90), 
          axis.title.x = element_blank())
}

pllist[[1]]
ggsave("data/aa34_svm_linear_conf_matrices_subst_grp_20190306.jpeg", height=7, width=7, units='in')
