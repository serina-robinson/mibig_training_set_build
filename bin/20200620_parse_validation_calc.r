# Install packages
pacman::p_load("tidyverse", "readxl", "DECIPHER")

x <- as.list(rnorm(10000))
names(x) <- paste("a", 1:length(x), sep = "")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the dataset
dat <- read_excel("output/Supplementary_Table_S3.xlsx") %>%
  janitor::clean_names()
dat
sum(dat$fc_correct)
sum(dat$ss_correct)

# Avg percent identity of fc
fc_sqs <- dat %>%
  dplyr::filter(fc_prediction_probability >= 0.6) %>%
  dplyr::summarise(summary(perc_aa_identity_to_closest_training_set_hit))
fc_sqs # perc aa 58.9


# Avg percent identity of ss
ss_sqs <- dat %>%
  dplyr::filter(ss_prediction_probability >= 0.6) %>%
  dplyr::summarise(summary(perc_aa_identity_to_closest_training_set_hit))
ss_sqs # perc aa 71.0, median 67
