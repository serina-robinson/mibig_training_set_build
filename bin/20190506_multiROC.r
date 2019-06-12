# Reauire packages
require(multiROC)
# install.packages('dummies')
require(dummies)
data(iris)
head(iris)

set.seed(123456)
total_number <- nrow(iris)
train_idx <- sample(total_number, round(total_number*0.6))
train_df <- iris[train_idx, ]
test_df <- iris[-train_idx, ]

rf_res <- randomForest::randomForest(Species~., data = train_df, ntree = 100)
rf_pred <- predict(rf_res, test_df, type = 'prob') 
rf_pred <- data.frame(rf_pred)
colnames(rf_pred) <- paste(colnames(rf_pred), "_pred_RF")
head(rf_pred)

mn_res <- nnet::multinom(Species ~., data = train_df)
mn_pred <- predict(mn_res, test_df, type = 'prob')
mn_pred <- data.frame(mn_pred)
head(mn_pred)
colnames(mn_pred) <- paste(colnames(mn_pred), "_pred_MN")

true_label <- dummies::dummy(test_df$Species, sep = ".")
true_label <- data.frame(true_label)
colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
colnames(true_label) <- paste(colnames(true_label), "_true")
final_df <- cbind(true_label, rf_pred, mn_pred)

final_df
roc_res <- multi_roc(final_df, force_diag=T)
pr_res <- multi_pr(final_df, force_diag=T)

plot_roc_df <- plot_roc_data(roc_res)
plot_pr_df <- plot_pr_data(pr_res)

require(ggplot2)
ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = Group, linetype=Method), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
               colour='grey', linetype = 'dotdash') +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))

ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = Group, linetype=Method), size=1.5) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.justification=c(1, 0), legend.position=c(.95, .05),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, 
                                         linetype="solid", colour ="black"))
