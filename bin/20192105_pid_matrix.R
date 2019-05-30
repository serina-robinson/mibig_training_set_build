## Install packages
pacman::p_load('Biostrings', 'DECIPHER')

## Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# # Read in the data
# dat <- readAAStringSet("data/669_training_set_to_extract_34aa_20192803.faa")
# 
# # Read in the sp2 A domains
# sp2dat <- readAAStringSet("data/sp2.adomains.faa")
# 
# # Combine the datasets
# comb <- c(dat, sp2dat)
# dedup <- comb[!duplicated(comb)]
# dedup2 <- dedup[!duplicated(names(dedup))]
# # 1759
# 
# # Write to file
# writeXStringSet(dedup2, "data/1759_combined_training_full_length_sqs.faa")


# Read in the pid matrix
pid <- fread("data/adenylpred_pid.txt", fill = T, data.table = F)
pid_num <- pid %>%
  select(-V1, -V2)
rowMeans(pid_num)
summary(na.omit(apply(pid_num, 1, min)))
which.min(na.omit(apply(pid_num, 1, min)))
na.omit(apply(pid_num, 1, min))

tmp <- na.omit(apply(pid_num, 1, min))
n <- length(tmp)
tmp
sort(tmp, partial = n - 1)[n-1]

min( tmp[tmp!=0.00] ) # 2.78% 
