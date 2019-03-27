pacman::p_load(
  "data.table",
  "ape",
  "ggtree",
  "ggplot2",
  "pamr",
  "RColorBrewer",
  "phangorn", 
  "ggthemes",
  "usethis",
  "Biostrings",
  "DECIPHER",
  "officer",
  "rvg",
  "drc",
  "tidyverse",
  "genoPlotR",
  "data.table",
  "janitor",
  "plot3D",
  "plot3Drgl",
  "scales",
  "readxl")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the adjusted alignment
aa <- readAAStringSet("data/713_aa34_signatures.fa")
names(aa)
phy<-read.phyDat("data/713_aa34_signatures.fa", format="fasta", type="AA")


# Make a distance matrix
# dm<-dist.ml(phy)
# mat <- as.matrix(dm)
# write.csv(mat, paste0("data/",numseqs,"_distance_matrix.csv"), quote = FALSE, row.names = T)
# write.table(mat, paste0("data/",numseqs,"_distance_matrix.tsv"), quote = FALSE, row.names = T, sep = "\t")

rd <- read_csv("data/713_seqs_136_feats_for_supervised.csv") 

rd2 <- rd %>%
  dplyr::select(-X1)
mds <- as.matrix(rd2)
rownames(mds) <- rd$X1
mdist <- dist(mds)

# Multidimensional scaling using cmdscale
mds<-cmdscale(mdist,eig=TRUE,k=3)
x<-mds$points[,1]
y<-mds$points[,2]
z<-mds$points[,3]

# Calculate percent explained by each principal component
pc1 <- mds$GOF[1]
pc2 <- (mds$GOF[2]-mds$GOF[1])
pc1
pc2

# Make x- and y-axis labels
xlab<-paste0("PC 1 (", round( pc1 * 100, 2), "% tot. explained var.)")
ylab<-paste0("PC 2 (", round( pc2 * 100, 2), "% tot. explained var.)")


# Plot 3D
# scatter3D(x,y,z)
# plotrgl()
# rgl::plot3d(x,y,z)

############# 
# Do sequences cluster by their substrate?
dat <- data.frame(cbind(x, y))
# rownames(dat) <- rd$X1
rownames(dat)
sub_tofind <- table(word(rownames(dat), -1, sep = "_"))


# Color points for MDS 
source("src/color_pcoa_points.r")
tocol <- sub_tofind
sub_tofind
pal <- palette(colorRampPalette(colors=brewer.pal(8,"Accent"))(length(sub_tofind)))

show_col(as.character(pal))

cdat <- color_pcoa_points(tocol, dat, pal)
head(cdat)
table(cdat$colrs)

# tr1 <- rownames(cdat)[cdat$y > 1 & cdat$nms == "long_chain"]
# tr2 <- rownames(cdat)[cdat$x < -0.5 & cdat$nms == "long_chain"]
# tr2
# tok <- str_split_fixed(tr1, "_", 3)[,c(1,2)]
# toko <- paste0(tok[,1], "_", tok[,2])
# tok2 <- str_split_fixed(tr2, "_", 3)[,c(1,2)]
# toko2 <-  paste0(tok2[,1], "_", tok2[,2])

# train2 <- read_excel("data/current_ANL_training_set.xlsx")
# train2$protein_names[train2$entry_name %in% toko]
# train2$protein_names[train2$entry_name %in% toko2]
# 
# 
# vlacs <- paste0(train2$entry_name[grep("very long|Very long", train2$protein_names)], "*")
# paste0(vlacs, collapse = "|")
# 
# vlac_pt <- cdat$nms[grep(paste0(vlacs, collapse = "|"), rownames(cdat))]
# vlac_lab <- grepl(paste0(vlacs, collapse = "|"), rownames(cdat))
# cdat$labl <- ifelse(vlac_lab, cdat$nms[vlac_lab], "")
# table(cdat$labl)

# Make a colored 2D plot
numseqs <- nrow(cdat)
numseqs
pdf(paste0("output/", numseqs,"_PCoA_colored.pdf"))
par(mar=c(0.01, 0.01, 0.01, 0.01))
ggplot(data = cdat, aes(x=x, y=y)) + 
  # geom_text(label = cdat$labl) +
  # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
  geom_point(fill = cdat$colrs, size = cdat$sz, 
             alpha = cdat$alph, shape = 21) +
  scale_color_manual(pal) +
  scale_shape(solid = TRUE) +
  labs(x=xlab,y=ylab) +
  theme_bw() +
  theme(axis.text = element_text(size = 14),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.position = c(0.14, 0.80),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title=element_text(size=20, face="bold", hjust=0)
  )
dev.off()

with(cdat, plot3d(x = x, y = y, z = z,
                  col = colrs, type = 's', size = 1))


# Try PCoA from the physicochemical properties?
