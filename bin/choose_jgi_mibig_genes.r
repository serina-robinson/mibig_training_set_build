## Install packages
pacman::p_load('ape', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set seed 
set.seed(123091)

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the MIBiG training data
rawdat <- read_excel("data/mibig_training_set_manually_edited_20190304.xlsx")

table(mibig$likely_substrate)
substrates_we_have <- c("(R)-3-hydroxymyristate", "salicylate", "butyrate", 
                        "propionate", "acetate", "cinnamate", "decanoate",
                        "hexanoate", "palmitate", "myristate", "octanoate",
                        "octanoate; hexanoate", "acetate; propionate", "acetate; butyrate; malonate",
                        "acetate; isobutyrate; propionate")
mibig <- rawdat %>%
  dplyr::filter(confidence < 4) %>%
  dplyr::filter(likely_substrate %in% substrates_we_have) %>%
  dplyr::filter(nchar(aa_seq) < 1000)
  # dplyr::filter(functional_class != "FAAL") %>%
  # dplyr::filter(!grepl("CP", transfer_molecule))
table(mibig$substrate_group)


write_csv(mibig, "candidate_mibig_JGI_genes.csv")
# only 34 with high confidence