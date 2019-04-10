## Install packages
pacman::p_load('ape', 'dendextend', 'readxl', 'rgl', 'cluster', 'ggrepel', 'plotly', 'tidyverse', 'data.table', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the aa classes
aaex <- read_excel("data/aa_cleaned_up_serina_min_grp_size_50.xlsx") %>%
  janitor::clean_names() %>%
 # dplyr::select(1:8) %>%
  map_df(tolower) %>%
  modify_at(1:8, ~ stringr::str_replace_all(., c("\\|" = "_")))
dim(aaex)

# map_df(~ gsub("-", "\\|",.x))
head(aaex)
colnames(aaex)

# Read in the adoms
adoms <- readAAStringSet("data/sp2.adomains.faa")
names(adoms) <- gsub("\\\t", "_", names(adoms))
names(adoms)
# adoms <- readAAStringSet("data/sp2_34extract_named.faa")
dtf <- data.frame(names(adoms), as.character(adoms), word(names(adoms), 2, sep = "_"), stringsAsFactors = F)
colnames(dtf) <- c("nms", "aa_seq", "subst")
head(dtf)

# Create a new column for larger amino acid class
dtf$subst <- tolower(dtf$subst)
table(dtf$subst)
dtf$subst <- gsub("\\|", "_", dtf$subst)
dtf$subst <- gsub("tyr_boh-tyr", "boh-tyr_tyr", dtf$subst)
dtf$subst <- gsub("tyr_phe", "phe_tyr", dtf$subst)
dtf$subst <- gsub("val_ile", "ile_val", dtf$subst)
dtf$subst <- gsub("^iso$", "ile", dtf$subst)
dtf$subst <- gsub("leu_arg", "arg_leu", dtf$subst)
dtf$subst <- gsub("pro_4me-pro", "4me-pro_pro", dtf$subst)
dtf$subst <- gsub("glu_3me-glu", "3me-glu_glu", dtf$subst)
dtf$subst <- gsub("hpg_2cl-hpg|hpg_cl2-hpg", "2cl-hpg_hpg", dtf$subst)
dtf$subst <- gsub("orn_oh-orn", "oh-orn_orn", dtf$subst)
dtf$subst <- gsub("trp_6cl-trp", "6cl-trp_trp", dtf$subst)
dtf$lrg_subst <- dtf$subst


for(j in 1:ncol(aaex))  {
  for(i in 1:nrow(dtf)) {
    if(dtf$subst[i] %in% pull(na.omit(aaex[,j]))) {
      dtf$lrg_subst[i] <- gsub(dtf$subst[i], colnames(aaex)[j], dtf$lrg_subst[i])
    }
  }
}

dtf$lrg_subst <- gsub("_", "\\.", dtf$lrg_subst)
head(dtf)
names(adoms) <- paste0(names(adoms), "_", dtf$lrg_subst)
head(names(adoms))
names(adoms) <- gsub(" ", "_", names(adoms))
names(adoms) <- gsub("____", "_", names(adoms))
table(word(names(adoms), sep ="_", -1))

tbdf <- data.frame(table(word(names(adoms), sep ="_", -1)))
tbdf <- tbdf[order(tbdf$Var1),]

all_members <- apply(aaex, 2, function(x) { paste0(na.omit(x), collapse = ", ") })
all_members <- gsub("_", "\\|", all_members)
all_members <- all_members[order(names(all_members))]
all_members
tbdf$members <- all_members
tbdf

write_csv(tbdf, "output/approx_large_aa_group_class_sizes.csv")
sub_dict <- data.frame(cbind(word(names(adoms), sep = "_", 2), dtf$lrg_subst), stringsAsFactors = F) %>%
  distinct()
dim(sub_dict)
head(sub_dict)
write_delim(col_names = F, sub_dict, "data/aa_substrate_class_dict.txt")

names(adoms)
# writeXStringSet(adoms, "data/1093_sp2_full_length_names_fixed_large_grps.faa")
# write_csv(data.frame(table(dtf$lrg_subst)), "output/not_included.csv")
# aa <- readAAStringSet("data/sp2_34extract_names_fixed_large_grps.faa")
# names(aa)


