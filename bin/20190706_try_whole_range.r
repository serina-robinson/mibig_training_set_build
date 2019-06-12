# Install packages
pacman::p_load("caret", "Biostrings", "phangorn", "ape", "seqinr", "DECIPHER", "cowplot", "tidymodels", "ranger", "tree", "rsample", "tidyverse", "randomForest","gbm","nnet","e1071","svmpath","lars","glmnet","svmpath")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the 34aa HMM aligned sequences
aa34 <- readAAStringSet("~/Documents/Wageningen_UR/github/check_a_dom_specificity/data/sp2_adoms_34extract_hmmalign.faa")
name_key <- readAAStringSet("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/data/sp2_34extract_names_fixed_large_grps.faa")
names(aa34) <- names(name_key)
table(word(names(aa34), sep = "_", 2))

# Remove all duplicates
dedup <- aa34[!duplicated(aa34)]
length(dedup) # 843


writeXStringSet(dedup, "data/843_34extract_hmmalign.faa")

# # Only take the first 9 residues
first9 <- AAStringSet(substr(dedup, 1, 9))
dedup <- first9[!duplicated(first9)]
length(dedup)
subst <- word(names(dedup), sep = "_", 2)
subst <- gsub("\\|", "\\.", subst)
subst <- gsub("\\-", "\\.", subst)
names(dedup) <- gsub("\\|", "\\.", names(dedup))
names(dedup) <- gsub("\\-", "\\.", names(dedup))

dtf <- data.frame(names(dedup), dedup, subst, stringsAsFactors = F) %>%
  dplyr::filter(!grepl("reject", names.dedup.)) %>%
  dplyr::filter(!grepl("-", dedup)) %>%
  group_by(subst) %>%
  add_count(subst) %>%
  dplyr::filter(n > 8)

dim(dtf) # 652
table(dtf$subst)
#
# # Remove all sequences with specificity less than 8
dedup2 <- AAStringSet(dtf$dedup)
# first9 <- AAStringSet(substr(dedup2, 1, 9))
table(duplicated(dedup2))
names(dedup2) <- dtf$names.dedup.
writeXStringSet(dedup2, "data/428_monomers_deduplicated.faa")

# Just train on the first 9 active site residues
rdaln <- read.alignment("data/428_monomers_deduplicated.faa", format = "fasta")

# Convert to a vector of physicochemical properties
source("src/convert_aln_15aap.r")
aa <- convert_aln_15aap(rdaln) #15 physicochemical properties
aadf <- data.frame(aa,stringsAsFactors = F)
dim(aadf)
write.csv(aadf, paste0("data/428_aa_seqs_135_feats_for_supervised_monomers_20190706.csv"), row.names=rownames(aadf), quote = F)

# Read in the data
rawdat <- read_csv("data/428_aa_seqs_135_feats_for_supervised_monomers_20190706.csv", skip_empty_rows = T)
rawdat <- data.frame(cbind(rawdat$X1), scale(rawdat[,2:ncol(rawdat)]), stringsAsFactors = F)
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, 2, sep = "_")
table(rawdat$clf)
dat <- rawdat[!duplicated(rawdat),]

rf_models <- list()
rf_pred <- list()
y_test_list <- list()
rf_cm <- list()

for(i in 1:50) {
  
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
  
  rf_models[[i]] <- ranger(y_train ~., data = form_train, num.trees = 1000, splitrule = "gini",
                           mtry = as.integer(sqrt(ncol(x_train))), min.node.size = 1,
                           importance = "permutation", class.weights = classwts)
  
  rf_pred[[i]] <- predict(rf_models[[i]], data = x_test)
  # rf_cm[[i]] <- confusionMatrix(rf_pred[[i]]$predictions, y_test)
  y_test_list[[i]] <- y_test
}


cm_res <- list()
for(i in 1:50) {
  tryCatch ({
    cm_res[[i]] <- confusionMatrix(rf_pred[[i]]$predictions, y_test_list[[i]])
  }, error=function(e){})
}
cm_res


# Calculate the overall test statistics
xx <- list()
ctr <- 1
cm_res[[1]]
for (i in 1:length(cm_res)) {
  if (is.null(cm_res[[1]])) {
    xx[[ctr]] <- cm_res[[i]]$overall[1]
    ctr <- ctr + 1
  }  
}

mean(unlist(xx))
unlist(xx)
sd(unlist(xx))
which.max(unlist(xx))

dtf <- data.frame(cm_res[[26]]$table)

# gsub("Udtf$Prediction
pdf("output/20190706_rf_monomers_cm.pdf", width = 8, height = 8)
ggplot(dtf, aes(Prediction, Reference)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90), 
        axis.title.x = element_blank(),
        legend.position = "none")
dev.off()

