## Install packages
pacman::p_load('ape', 'tidyverse', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set seed 
set.seed(123091)

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the MIBiG training data
rawdat <- read_excel("data/mibig_training_set_manually_edited_20192603.xlsx") 
rawdat$large_substrate_group <- gsub("_", "", rawdat$large_substrate_group)
rawdat$functional_class <- gsub("PEPTIDE", "NRPS", rawdat$functional_class)
rawdat$small_substrate_group[rawdat$likely_substrate == "phenylacetate"] <- "cinnamate.and.succinylbenzoate.derivatives"
write_csv(rawdat, "data/mibig_training_set_manually_edited_20192603.csv")

mibig <- rawdat %>%
  # dplyr::filter(confidence > 0) %>%
  mutate(substrate_group = str_replace_all(substrate_group, "coumarin", "aryl")) %>%
  mutate(substrate_group = str_replace_all(substrate_group, "biaryl", "aryl")) %>%
  mutate(substrate_group_tr = str_replace_all(substrate_group, "_", "")) %>%
  mutate(functional_class = str_replace_all(functional_class, "BIARYL", "ARYL")) %>%
  mutate(functional_class = str_replace_all(functional_class, "COUM", "ARYL")) %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  mutate(sqnams_tr = paste0(bgcs, "_", word(cmpnd, 1, sep = "_"), "_", acc, "_", substrate_group_tr, "_", functional_class)) %>%
  mutate(sqnams_tr = str_replace_all(sqnams_tr, "-", "_")) %>%
  mutate(sqnams_tr = str_replace_all(sqnams_tr, "\\.", "_"))

mib2 <- mibig %>%
  dplyr::rename(large_substrate_group = substrate_group) %>%
  mutate(small_substrate_group = case_when(functional_class == "OTHER|HOLDOUTTEST" ~ "unknown.other",
                                #  str_detect(likely_substrate, "phenylacet") ~ "phenylacetate",
                                  str_detect(likely_substrate, "carnitine|isoleucine|leuc|proli|alanine|aminolev|threonine") ~ "amino.acid",
                                  str_detect(likely_substrate, "stear|olei|olea|docosa|tetracosa|oxocholes|nonadec|arachi|eicosa|bile|verylongchain") ~ "C18.and.up.or.bile.acid",
                                  # str_detect(likely_substrate, "pyrrol|pyrol") ~ "pyrrole_derivative",
                                  str_detect(likely_substrate, "tridec|myrist|palmit|longchain|phytodien|heptadec|long_chain") ~ "C13.through.C17",
                                  str_detect(likely_substrate, "pent|valer|hex|hept|oct|mediumchain") ~ "C5.through.C8",
                                  str_detect(likely_substrate, "non|dec|dodec|laur") ~ "C9.through.C12",
                                  str_detect(likely_substrate, "acetat|shortchain|acryl|prop|but|crot") ~ "C2.through.C4",
                                  # str_detect(likely_substrate, "chlorobenz") ~ "halogenated_benzyl_derivatives",
                                  str_detect(likely_substrate, "dihydroxybenz|salicyl|picolin|nicotin|phenylacet") ~ "salicylic.acid.derivatives",
                                  str_detect(likely_substrate, "coum|cinna|vanill|ferul|succinyl|quinol") ~ "cinnamate.and.succinylbenzoate.derivatives",
                                  str_detect(likely_substrate, "luciferin") ~ "luciferin",
                                  str_detect(likely_substrate, "medianbetahydroxyacidlong") ~ "median.beta.hydroxyacid",
                                  str_detect(likely_substrate, "anthran|benzo|naphth|aryl|xanthu|quin") ~ "aryl.and.biaryl.derivatives",
                                  str_detect(likely_substrate, "malon|methylmal|oxala|acetoace") ~ "malonate.derivatives",
                                  TRUE ~ "unknown.other"))
# write_csv(mib2, "output/mibig_training_set_manually_edited_20192603.csv")

mib3 <- read_excel("data/mibig_training_set_manual_edited_20192603.xlsx")

# Read in the UniPROT KB training data
# uniprot <- read_csv("data/anl_training_set_updated_20192503_fixnams.csv") %>%
#   dplyr::filter(!grepl("coelente", substrate)) %>%
#   mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
#   mutate(sqnams_tr = paste0(1:nrow(.), "_", org_short, "_", substrate_group, "_", functional_class))
# uniprot$sqnams_tr

# Combine the MIBiG and the UNIPROT
# dat <- uniprot %>%
#   dplyr::mutate(likely_substrate = substrate) %>%
#   dplyr::mutate(pmid = pub_med_id) %>%
#   bind_rows(., mibig) %>%
#   dplyr::select(entry_name, organism, protein_names, likely_substrate, substrate_group, functional_class,
#                 bgcs, pmid, cmpnd, pdb_id, kinetics, ec_numbers, title, sqnams_tr, aa_seq) #%>%
# 


# Reclassify likely substrate
datsub <- dat %>%
  # dplyr::filter(functional_class != "HOLDOUTTEST") %>%
  mutate(small_subclass = case_when(functional_class == "OTHER|HOLDOUTTEST" ~ "unknown-other",
                                    str_detect(likely_substrate, "phenylacet") ~ "phenylacetate",
                                    str_detect(likely_substrate, "carnitine|isoleucine|leuc|proli|alanine|aminolev|threonine") ~ "amino-acid",
                                    str_detect(likely_substrate, "stear|olei|olea|docosa|tetracosa|oxocholes|nonadec|arachi|eicosa|bile|verylongchain") ~ "C18+-or-bile-acid",
                                    # str_detect(likely_substrate, "pyrrol|pyrol") ~ "pyrrole_derivative",
                                    str_detect(likely_substrate, "tridec|myrist|palmit|longchain|phytodien|heptadec|long_chain") ~ "C13-C17",
                                    str_detect(likely_substrate, "pent|valer|hex|hept|oct|mediumchain") ~ "C5-C8",
                                    str_detect(likely_substrate, "non|dec|dodec|laur") ~ "C9-C12",
                                    str_detect(likely_substrate, "acetat|shortchain|acryl|prop|but|crot") ~ "C2-C4",
                                    # str_detect(likely_substrate, "chlorobenz") ~ "halogenated_benzyl_derivatives",
                                    str_detect(likely_substrate, "dihydroxybenz|salicyl|picolin|nicotin") ~ "salicylic-acid-derivatives",
                                    str_detect(likely_substrate, "coum|cinna|vanill|ferul|succinyl|quinol") ~ "cinnamate-and-succinylbenzoate-derivatives",
                                    str_detect(likely_substrate, "luciferin") ~ "luciferin",
                                    str_detect(likely_substrate, "medianbetahydroxyacidlong") ~ "median-beta-hydroxyacid",
                                    str_detect(likely_substrate, "anthran|benzo|naphth|aryl|xanthu|quin") ~ "aryl-and-biaryl-derivatives",
                                    str_detect(likely_substrate, "malon|methylmal|oxala|acetoace") ~ "malonate-derivatives",
                                    TRUE ~ "unknown-other"))


table(datsub$small_subclass)
datseq <- AAStringSet(datsub$aa_seq)
names(datseq) <- paste0(datsub$small_subclass, "_", datsub$sqnams_tr)
names(datseq)

# dtf <- data.frame(cbind(datsub$likely_substrate, datsub$small_subclass), stringsAsFactors = F)
#write_csv(dtf, "output/substrates_and_likely_substrate_classes.csv")

nrps <- readAAStringSet("../sandpuma2_serina/flat/sp2.adomains.faa")
monomers <- word(names(nrps), 2, sep = "\\\t")
monomers
ids <- gsub("\\.", "_", word(names(nrps), 3, sep = "\t"))
ids
names(nrps) <- paste0("amino-acid","_",ids, "_", monomers, "_aminoacid_NRPS")
# names(nrps) <- paste0(monomers, "_", names(nrps), "_aminoacid_NRPS")
names(nrps) <- gsub("\\\t", "_", names(nrps))
names(nrps)

mondf <- data.frame(cbind(monomers, names(nrps)), stringsAsFactors = F) %>%
  dplyr::group_by(monomers) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()
mondf$monomers
table(mondf$monomers) # take one from each monomer
nrps_rep <- AAStringSet(nrps[names(nrps) %in% mondf$V2])
nrps_rep

sqs <- c(datseq, nrps_rep)
head(names(sqs))

# Remove duplicates
nodups <- sqs[!duplicated(sqs)]
dups <- sqs[duplicated(sqs)]
names(dups)
length(nodups) # 3 duplicates
nodups

# Remove sequences shorter than 400 amino acids
sqs_long <- nodups[width(nodups) >= 350] # 575 sequences
sqs_short <- nodups[width(nodups) < 350]
length(sqs_long) # 740 left
summary(width(sqs_long))
length(sqs_long)
table(word(names(sqs_long), sep = "_", 1))
writeXStringSet(sqs_long, "../sandpuma2_serina/flat/740_training_seqs_small_subclass.fa")

