## Install packages
pacman::p_load("caret", "Biostrings", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
rawdat <- readAAStringSet("../sandpuma2_serina/flat/740_training_set_34extract.faa")
length(rawdat)
names(rawdat)

# Convert the 740 aa signatures to features
rdaln <- read.alignment(file = '../sandpuma2_serina/flat/740_training_set_34extract.faa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)
aa <- bgafun::convert_aln_AAP(rdaln) #5 physicochemical properties

aadf <- data.frame(aa, stringsAsFactors = F)

aap <- aadf %>%
  dplyr::mutate(nms = rownames(.)) %>%
  dplyr::select(-contains("D")) %>%
  dplyr::filter(!grepl("ERROR", X1A))

rownames(aap) <- aap$nms
aap <- aap %>%
  dplyr::select(-nms)
colnames(aap) <- gsub("^X","",colnames(aap))

colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap))
numfeats <- length(colnames(aap)) # 1096
colnames(aap)

# Write to CSV file
# write.csv(aap, paste0("data/740_seqs_", numfeats, "_feats_for_supervised.csv"), row.names=rownames(aap), quote = F)

# Read in the data
rawdat <- read_csv(paste0("data/740_seqs_", numfeats, "_feats_for_supervised.csv"))
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 1, sep = "_")
rawdat$clf <- gsub("-", "\\.", rawdat$clf)
rawdat$clf <- gsub("\\+", ".and.up", rawdat$clf)
table(rawdat$clf)

# Remove the holdout test predictions
dat <- rawdat[-grep(paste0(c("unknown-other"), collapse = "|"), rawdat$nms),] # 658 observations

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
  preProc = c("center", "scale"),
  trControl = trainControl(method = "cv", number = 10,
                           verboseIter = T, classProbs = T,
                           savePredictions = "final"),
  # tuneGrid = tgrid,
  num.trees = 500,
  # respect.unordered.factors = FALSE,
  verbose = TRUE,
  importance = "permutation")

# 
warnings()
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
# form_train <- data.frame(cbind(form_train[,1],
#                                scale(form_train[,2:(ncol(form_train)-1)], center = T),
#                                form_train[,ncol(form_train)]), stringsAsFactors = F)
# rf_ml <- ranger(y_train ~., data = form_train, num.trees = 500, splitrule = "gini", mtry = 69, importance = "permutation")
mod <- readRDS("output/rf_small_subst_grp_20192603.rds")
# makpred <- predict(mod, newdata = aap)

rf_pred <- predict(rf, newdata = form_test)
cm_rf <- confusionMatrix(rf_pred, as.factor(dat_test$clf))
sink("data/cm_rf_substrate_small_grp.txt")
cm_rf
sink()


## Make a heatmap of confusion matrix results

cm_list <- list(cm_rf$table)
names(cm_list) <- c("rf_aa_34")
colnames(data.frame(cm_list[[1]]))

str(cm_list[[1]])
pllist <- list()
for(i in 1:length(cm_list)) {
  pllist[[i]] <- ggplot(data.frame(cm_list[[1]]), aes(Prediction, Reference)) +
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

ggsave("aa34_group_conf_matrices.jpeg", height=8, width=8, units='in')
# saveRDS(rf, "output/rf_small_subst_grp_20192603.rds")
