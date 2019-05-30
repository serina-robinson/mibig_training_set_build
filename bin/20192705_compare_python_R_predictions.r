# Install packages
pacman::p_load("data.table", "tidyverse", "plotKML", "scales", "tidyverse", "stringr", "Biostrings","DECIPHER","igraph","RColorBrewer", "officer", "rvg")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the Python predictions
colnames(pypred)
pypred <- fread("data/parsed_BGC_sequences_combined_27052019_predictions.txt", data.table = F) %>%
  janitor::clean_names() %>%
  mutate(nam_trimmed = substr(query_name, 1, (nchar(query_name) - 5))) %>%
  dplyr::select(nam_trimmed, prediction, probability_score) %>%
  dplyr::rename(query_name = nam_trimmed,
                py_prediction = prediction,
                py_probability_score = probability_score)

pypred$py_prediction <- gsub("\\.", " ", pypred$py_prediction)
pypred$py_prediction <- gsub(" aa", "", pypred$py_prediction)

# Read in the R predictions
rpred <- read_csv("data/AdenylPred-2019-05-27.csv") %>%
  janitor::clean_names()
rpred$predicted_substrate_specificity_ss <- tolower(rpred$predicted_substrate_specificity_ss)


# Join the two predictions
comb_pred <- rpred %>%
  inner_join(., pypred, by = "query_name")
head(comb_pred)
write_csv(comb_pred, "output/combined_py_and_r_predictions.csv")

# See which agree
agree_pred <- comb_pred[comb_pred$py_prediction == comb_pred$predicted_substrate_specificity_ss,] # most agree
dont_agree <- comb_pred[comb_pred$py_prediction != comb_pred$predicted_substrate_specificity_ss,]
write_csv(dont_agree, "output/dont_agree_combined_py_and_r_predictions.csv")

