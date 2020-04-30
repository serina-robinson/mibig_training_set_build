# Install packages
pacman::p_load("webchem", "readxl", "tidyverse", "ChemmineR", "data.table", "plot3D", "plot3Drgl",
               "ggrepel")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("ChemmineR", version = "3.8")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in all the chemical structures from the Uniprot training set
uniprot <- read_excel('data/anl_training_set_updated_20190215.xlsx') %>%
  dplyr::filter(!grepl("coelente", substrate)) %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  mutate(org_short = gsub(" ", "_", paste0(word(organism, sep = " ", 1), "_", word(organism, sep = " ", 2)))) %>%
  mutate(sqnams_tr = paste0(1:nrow(.), " ", org_short, "_", substrate_group, "_", functional_class))

# Read in the MIBiG training data
rawdat <- read_excel("data/mibig_training_set_manually_edited_20190304.xlsx")
mibig <- rawdat %>%
  dplyr::filter(confidence > 0) %>%
  mutate(substrate_group = str_replace_all(substrate_group, "coumarin", "aryl")) %>%
  mutate(substrate_group = str_replace_all(substrate_group, "biaryl", "aryl")) %>%
  mutate(substrate_group_tr = str_replace_all(substrate_group, "_", "")) %>%
  mutate(functional_class = str_replace_all(functional_class, "BIARYL", "ARYL")) %>%
  mutate(functional_class = str_replace_all(functional_class, "COUM", "ARYL")) %>%
  mutate(functional_class = str_replace_all(functional_class, "MMCS", "SACS")) %>%
  mutate(sqnams_tr = paste0(bgcs, "_", word(cmpnd, 1, sep = "_"), "_", acc, "_", substrate_group_tr, "_", functional_class)) %>%
  mutate(sqnams_tr = str_replace_all(sqnams_tr, "-", "_")) %>%
  mutate(sqnams_tr = str_replace_all(sqnams_tr, "\\.", "_")) %>%
  mutate(substrate = likely_substrate)

# Combine uniprot and mibig

sub_grps <- c("very_long_chain", "long_chain", "median_beta_hydroxy_acid_long", "aryl_polyene", "aryl",
              "medium_chain", "unknown", "C24_bile_acids", "C27_bile_acid", "bile_acid", "short_chain") 

cmbnd <-  uniprot %>%
  bind_rows(mibig, .) %>%
  dplyr::filter(!grepl(paste0(sub_grps, collapse = "|"), substrate)) %>%
  mutate(spec_substrate = gsub("_", " ", substrate))
# dplyr::select(spec_substrate)   
sub_key <- data.frame(cmbnd$likely_substrate, cmbnd$substrate_group, stringsAsFactors = F)
head(sub_key)

# sum(is.na(cmbnd$substrate)) # no NAs
spec_substrate <- gsub("\\/", ";", cmbnd$spec_substrate)
spec_substrate <- gsub("nicotinic acid", "nicotinoate", spec_substrate)
res <- str_split(pattern = ";", string = trimws(spec_substrate), simplify = TRUE)
res_df <- data.frame(trimws(res))
rownames(res_df) <- cmbnd$sqnams_tr
tail(res_df)
# spec_substrate_unl <- unlist(lapply(1:length(spec_substrate), function(x) { str_split(";", trimws(spec_substrate[x])) }))
# write_csv(res_df, "data/proteins_specific_substrates_grouped.csv", row.names = T)
table(res_df$X2)

test_pull <- unique(unlist(sapply(1:ncol(res_df), function(x) as.character(res_df[,x])[res_df[,x] != ""])))
# all_cids <- get_cid(test_pull)


cid_df <- data.frame(as.matrix(all_cids)) %>%
  rename(cids = as.matrix.all_cids.) %>%
  unlist() %>%
  na.omit() %>%
  as.data.frame()

cid_df$cmpnd <- rownames(cid_df)
final_df <- cid_df %>%
  dplyr::filter(!grepl("2$", cmpnd)) 

colnames(final_df)
colnames(final_df)[1] <- "cid"
final_df <- final_df %>%
  mutate(cid = as.numeric(cid))

final_df
final_df$cid
# write_delim(data.frame(final_df$cid), "data/all_cids_output.txt", delim = "\t", col_names = F)
final_df <- read_delim("data/all_cids_output.txt", delim = "\t", col_names = F)
head(final_df)
# https://pubchem.ncbi.nlm.nih.gov/pc_fetch/pc_fetch.cgi

# Read in all the chemical structures from the MIBiG training set 
smiles <- fread(file = "data/cids_SMILES_key.txt", sep = "\t", header = F, data.table = F)

# Read in the SMILES 
# Read in the SDF structures
sdfset <- read.SDFset("data/90_cids_sdf.sdf")
cid_pulled <- as.numeric(unlist(lapply(1:length(sdfset), function(x) { header(sdfset[[x]])[1] })))

cid_joined <- cid_pulled %>%
  as.data.frame(as.integer(.))

colnames(cid_joined) <- "cid"

colnames(final_df) <- "cid"
final_join <- final_df %>%
  dplyr::filter(!duplicated(cid)) %>%
  inner_join(., cid_joined, by = "cid", all.x = FALSE) 
final_join
# write_csv(final_join, "data/final_join_cmpnd_cid_key.csv")

final_join <- read_csv("data/final_join_cmpnd_cid_key.csv")

apset <- sdf2ap(sdfset)
apset@ID <- gsub("cids.", "", final_join$cmpnd)
apset@ID <- gsub("1$", "", apset@ID)

# keys <- read_delim("data/all_cids_output.txt", delim = "\t")

# apdups <- cmp.duplicated(apset, type=1)
# # sdfset2 <- sdfset[which(!apdups)]
# # apset2 <- apset[which(!apdups)]
# # fpset2 <- desc2fp(apset2)


clusters <- cmp.cluster(db=apset,  cutoff=0.41, save.distances = "data/distmat.rda")
head(clusters)
write_csv(clusters, "output/tanimoto_clusters.csv")
cluster.sizestat(clusters) # 56, 25 ,2


clusviz <- cluster.visualize(apset, clusters, size.cutoff = 1, quiet = TRUE)
x <- clusviz[,1]
y <- clusviz[,2]
trdat <- data.frame(cbind(x,y))
rownames(trdat) <- names(x)
trdat$likely_substrate <- rownames(trdat)
colnames(sub_key) <- gsub("cmbnd\\.", "", colnames(sub_key))

head(trdat)
dim(trdat)

# dim(trdat[trdat$likely_substrate %in% sub_key$likely_substrate,])

subgrp <- trdat %>%
  left_join(., sub_key, by = "likely_substrate") %>%
  distinct()
# distinct()

# write_csv(subgrp, "output/likely_substrate_for_mds.csv")

subgrp <- read_csv("output/likely_substrate_for_mds.csv")
is.na(subgrp$substrate_group)

dev.off()
subgrp$likely_substrate[grepl("carboline|naphthoate", subgrp$likely_substrate)] <- ""
pdf(paste0("output/", numseqs,"_PCoA_labeled.pdf"), width = 12, height = 12)
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = subgrp, aes(x=x, y=y, label = likely_substrate)) + 
  # geom_text(label = trdat$nms) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  geom_point() +
  geom_point(aes(fill = as.factor(subgrp$substrate_group)), shape = 21, size = 4) +
  scale_shape(solid = TRUE) +
  geom_text_repel() +
  # labs(x=xlab,y=ylab) +
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = "top",
        legend.box = "horizontal",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=14, face="bold", hjust=0))
dev.off()

library(scatterplot3d) 
coord <- cluster.visualize(apset, clusters, size.cutoff=1, dimensions=3, quiet=TRUE) 
library(rgl) 
rgl.open(); offset <- 50;
par3d(windowRect=c(offset, offset, 640+offset, 640+offset)) 
rm(offset)
rgl.clear() 
rgl.viewpoint(theta=45, phi=30, fov=60, zoom=1)
spheres3d(coord[,1], coord[,2], coord[,3], radius=0.03, color=coord[,4], alpha=1, shininess=20) 
aspect3d(1, 1, 1) 
axes3d(col='black')
title3d("", "", "", "", "", col='black')
bg3d("white") 

load("data/distmat.rda") 
coord <- cluster.visualize(apset, clusters, size.cutoff=1, dimensions=3, quiet=TRUE) 


p <- plot_ly(cdat, x = Comp.1 , y = Comp.2, text = rownames(carsDf),
             mode = "markers", color = cluster_name, marker = list(size = 11))