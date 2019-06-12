# Install packages
pacman::p_load("tidyverse")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the data
dat <- read_csv("data/combined_adenylate_forming_training_set_20192903.csv")
colnames(dat)

coum <- dat %>%
  dplyr::filter(grepl("coum|cinna|vanill|ferul", substrate)) %>%
  dplyr::select(-likely_substrate)

write_csv(coum, "output/Coumarate_CoA_ligases.csv")
