# Install packages
pacman::p_load("DECIPHER", "tidymodels", "tidyverse", "readxl")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Set seed 
set.seed(20190304)

# Read in the data
rawdat <- read_csv("data/1553_training_sqs_with_loop_extracted.csv")
colnames(rawdat)[1] <- "nms"
rawdat$clf <- word(rawdat$nms, -1, sep = "_")
table(rawdat$clf)
rawdat$nms[rawdat$clf == "luciferin"]


# Find HOLDOUT sequences
init_dat <- read_csv("../antismash_clustering/data/6052_seqs_1096_feats_for_supervised.csv")
holdouts <- init_dat[grep("HOLDOUT", init_dat$X1),]

# Remove those for which the true substrate is not known
holdouts_fin <- holdouts %>%
  dplyr::filter(!grepl("unknown|medianbetahydroxyacidlong", X1))
table(word(holdouts_fin$X1, -2, sep = "_")) # 8 aryl, # 7 long chain, 10 short chain
# Need luciferases, very long chain, if possible, and NRPS
bgcs <- word(holdouts_fin$X1, sep = "_", 3)
accs <- word(holdouts_fin$X1, sep = "_", 3)
accs[1] <- word(holdouts_fin$X1[1], sep = "_", 5)

# Check if holdout predictions are in training set
holdcheck <- rawdat$nms[grep(paste0(accs, collapse = "|"), rawdat$nms)]
holdcheck
holdcheck_accs <- word(holdcheck, 1, sep = "\\.1")
leftout <- holdouts_fin$X1[!grepl(paste0(holdcheck_accs, collapse = "|"), holdouts_fin$X1)]
table(word(leftout, sep = "_", -2))
leftout_accs <- word(leftout, sep = "_", -4)

# Search for leftout in MIBIG entry
mibig <- read_excel("data/mibig_training_set_manually_edited_20192603.xlsx")
leftout_accs
which_mib <- mibig[grep(paste0(leftout_accs, collapse = "|"), mibig$acc),]
dim(which_mib)
mibaa <- AAStringSet(which_mib$aa_seq)
names(mibaa) <- leftout
names(mibaa)
leftout
# writeXStringSet(mibaa, "output/21_holdout_seqs_to_test_20200414.fasta")


# Include others
dat_other <- read_excel("data/combined_adenylate_forming_training_set_for_db_20191404.xlsx") %>%
  dplyr::filter(functional_class == "OTHER")
other_accs <- word(dat_other$acc, sep = "\\.1", 1)
found <- rawdat$nms[grep(paste0(other_accs, collapse = "|"), rawdat$nms)]
otheraa <- AAStringSet(dat_other$aa_seq)
names(otheraa) <- paste0(other_accs, "_", dat_other$organism, "_", dat_other$substrate)
names(otheraa) <- gsub("\\(|\\)|\\[|\\]| |\\/", "_", names(otheraa))

comball <- AAStringSet(c(otheraa, mibaa))

reffind <- 
# writeXStringSet(comball,"output/35_weird_ANL_members.fasta")


# Pull new MIBiG entries?
# sp2 <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
# bgcs <- word(names(sp2), sep = "BGC", 2)
# bgc2 <- word(bgcs, sep = "_", 1)
# bgc3 <- paste0("BGC", bgc2[!is.na(bgc2)])
# bgcsort <- sort(bgc3, decreasing = T)
# bgcsort # Last BGC included was 1479

# tab1 <- read_excel("data/combined_adenylate_forming_training_set_for_db_20191404.xlsx")
# tab2 <- read_excel("data/combined_adenylate_forming_training_set_20192803.xlsx")
# tab3 <- read_csv("data/combined_adenylate_forming_training_set_20192703.csv")
