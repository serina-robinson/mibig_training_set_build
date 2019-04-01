## Install packages
pacman::p_load('ape', 'dendextend', 'rgl', 'cluster', 'ggrepel', 'plotly', 'tidyverse', 'data.table', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the training set
adoms <- readAAStringSet("data/sp2_34extract_named.faa")
names(adoms)

# A-domain classes 
subst <- word(names(adoms), sep = "_", 2)
subst_nms <- paste0(names(adoms), "_NRPS")
length(table(word(subst_nms, 2, sep = "_"))) # 123 classes

# PCA of all amino acids
# Convert the 713 aa signatures to features
rdaln <- read.alignment(file = 'data/sp2_34extract_named.faa', format = "fasta")
rdaln$seq <- toupper(rdaln$seq)
aa <- bgafun::convert_aln_AAP(rdaln) #5 physicochemical properties

aadf <- data.frame(aa, stringsAsFactors = F)

aap <- aadf %>%
  dplyr::mutate(nms = rownames(.)) %>%
  dplyr::select(-contains("D")) %>%
  dplyr::filter(!grepl("ERROR", X1A))

rownames(aap) <- aap$nms
aap <- aap %>%
  dplyr::select(-nms)
colnames(aap) <- gsub("^X","",colnames(aap))

colnames(aap) <- paste0(c("polrty", "secstr", "molsz", "elechrg"), "_", colnames(aap))
numfeats <- length(colnames(aap)) # 136 features
colnames(aap)
numfeats


# write_csv(aap, "data/1093_aa_136_feats.csv")
rd <- read_csv("data/1093_aa_136_feats.csv")
rdscale <- scale(rd, center = TRUE)
head(rdscale)
# Sequence based dist
# phy<-read.phyDat("data/sp2_34extract_named.faa", format="fasta", type="AA")

#Compute a distance matrix
# sdist<-dist.ml(phy)

# Physicochemical dist
mds <- as.matrix(rd)
rownames(mds) <- rownames(aap)
mdist <- dist(mds)

# Multidimensional scaling using cmdscale
mds<-cmdscale(mdist, eig=TRUE, k=136)
x<-mds$points[,1]
y<-mds$points[,2]
z<-mds$points[,3]



# Calculate percent explained by each principal component
pc1 <- mds$GOF[1]
pc2 <- (mds$GOF[2]-mds$GOF[1])
pc1
pc2

# Make x- and y-axis labels


# Plot 3D
# scatter3D(x,y,z)
# plotrgl()
# rgl::plot3d(x,y,z)

############# 
# Do sequences cluster by their substrate?
dtf <- data.frame(names(adoms), as.character(adoms), word(names(adoms), 2, sep = "_"), stringsAsFactors = F)
colnames(dtf) <- c("nms", "aa_seq", "subst")

nrow(dtf)
dat <- data.frame(cbind(x, y, dtf$subst), stringsAsFactors = F)
# rownames(dat) <- rd$X1
rownames(dat)
sub_tofind <- table(word(rownames(dat), 2, sep = "_"))
sub_tofind

# 10 or more
ten_or_more <- dtf %>%
  add_count(subst) %>%
  dplyr::filter(n > 9)
# ten_or_more[order(ten_or_more$n, decreasing = T),]


# Color points for MDS 
source("src/color_pcoa_points.r")
tocol <- unique(ten_or_more$subst)
length(sub_tofind)
pal <- palette(colorRampPalette(colors=brewer.pal(8,"Accent"))(length(sub_tofind)))

show_col(as.character(pal))
cdat <- color_pcoa_points(tocol, dat, pal)
table(cdat$colrs)
  
# Make a colored 2D plot
numseqs <- nrow(cdat)
cdat$x <- as.numeric(cdat$x) 
cdat$y <- as.numeric(cdat$y)
cdat$labels <- rep("", nrow(cdat))
cdat$labels
for(i in 1:length(unique(cdat$V3))) {
  inds <- grep(unique(cdat$V3)[i], cdat$V3)[1]
  cdat$labels[inds] <- unique(cdat$V3)[i]
}


pl <- list()
for(i in 1:5) {
  pl[[i]] <- list()
  for(j in 1:5) {
  cdat$x <- mds$points[,i]
  cdat$y <- mds$points[,j]
  xlab<-paste0("PC ", i, " (", round( pc1 * 100, 2), "% tot. explained var.)")
  ylab<-paste0("PC ", j, " (", round( pc2 * 100, 2), "% tot. explained var.)")
  # pdf(paste0("output/", numseqs,"_PCoA_colored.pdf"), width = 10, height = 10)
  # par(mar=c(0.01, 0.01, 0.01, 0.01))
  pl[[i]][[j]] <- ggplot(data = cdat, aes(x=x, y=y)) + 
    geom_text_repel(label = cdat$labels) +
    # geom_text(label = cdat$labl) +
    # geom_text(x = cdat$x, y = cdat$y, label = cdat$labl, check_overlap = T) +
    geom_point(fill = pal[as.numeric(as.factor(cdat$V3))], size = 3, shape = 21) + 
               # alpha = cdat$alph, shape = 21) +
    scale_shape(solid = TRUE) +
    labs(x=xlab,y=ylab) +
    theme_bw() +
    theme(# axis.text = element_text(size = 14),
          # legend.key = element_rect(fill = "white"),
          # legend.background = element_rect(fill = "white"),
          # legend.position = c(0.14, 0.80),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # plot.title=element_text(size=20, face="bold", hjust=0)
    )
  }
}




# pdf("data/PCoA_aa_plot.pdf", width = 50, height = 50)
plot_grid(pl[[1]][[1]], pl[[1]][[2]], pl[[1]][[3]], pl[[1]][[4]], pl[[1]][[5]],
          pl[[2]][[1]], pl[[2]][[2]], pl[[2]][[3]], pl[[2]][[4]], pl[[2]][[5]],
          pl[[3]][[1]], pl[[3]][[2]], pl[[3]][[3]], pl[[3]][[4]], pl[[3]][[5]],
          pl[[4]][[1]], pl[[4]][[2]], pl[[4]][[3]], pl[[4]][[4]], pl[[4]][[5]],
          pl[[5]][[1]], pl[[5]][[2]], pl[[5]][[3]], pl[[5]][[4]], pl[[5]][[5]], 
          nrow = 5, ncol = 5)
# dev.off()


with(cdat, plot3d(x = x, y = y, z = z,
                  col = colrs, type = 's', size = 1))

# Currently - extract only one per class

# 
# nrow(dtf) - nrow(ten_or_more) # 194 that trimmed

# p <- plot_ly(data = cdat, x = x , y = y, showlegend = FALSE,
#              mode = "markers", color = cdat$colrs, marker = list(size = 11),
#              text = rownames(cdat), hoverinfo = 'text')
# p
# Find these in the sp2
# sp2 <- adoms[names(adoms) %in% onedf$X2]
# writeXStringSet(sp2, "data/adom_signatures_one_per_class_34aa_20190328.faa")

# Convert
# Make a legend
# pdf(file=paste0("output/",nrow(net),"_legend_aa34.pdf"))
# plot.new()
# legend("bottomright", legend=proteins, fill = colors, bty="n")
# dev.off()

# PAM medoids (Partitioning around medoids)

pm <- pam(mat, k = 10, metric = "euclidean", stand = FALSE)

dtf10 <- data.frame(rownames(cdat), names(pm$clustering),  pm$clustering)
colnames(dtf10) <- c("full_id", "substrate", "cluster_id")
head(dtf10)
# write_csv(dtf10, "output/amino_acid_10_clustering.csv")

pm <- pam(mat, k = 15, metric = "euclidean", stand = FALSE)
dtf15 <- data.frame(rownames(cdat), names(pm$clustering),  pm$clustering)
colnames(dtf15) <- c("full_id", "substrate", "cluster_id")
# write_csv(dtf15, "output/amino_acid_15_clustering.csv")

pm <- pam(mat, k = 20, metric = "euclidean", stand = FALSE)
dtf20 <- data.frame(rownames(cdat), names(pm$clustering),  pm$clustering)
colnames(dtf20) <- c("full_id", "substrate", "cluster_id")
table(dtf15$cluster_id)
# write_csv(dtf20, "output/amino_acid_20_clustering.csv")

# PAM medoids (Partitioning around medoids)
# pm <- pam(mdist, k = 8, metric = "euclidean", stand = FALSE)
# dtf8 <- data.frame(names(pm$clustering), cdat$V3, pm$clustering)
# colnames(dtf8) <- c("full_id", "substrate", "cluster_id")
# write_csv(dtf8, "output/amino_acid_8_clustering.csv")
# 
# # PAM medoids (Partitioning around medoids)
# pm <- pam(mdist, k = 9, metric = "euclidean", stand = FALSE)
# dtf9 <- data.frame(names(pm$clustering), cdat$V3, pm$clustering)
# colnames(dtf9) <- c("full_id", "substrate", "cluster_id")
# write_csv(dtf9, "output/amino_acid_9_clustering.csv")

dtftabs <- dtf10 %>%
  group_by(cluster_id, substrate) %>%
  add_count(substrate) 
dtftabs  
  
clusterdf <- data.frame(table(dtf10$cluster_id, dtf10$substrate))
clusterfilt <- clusterdf %>%
  dplyr::filter(Freq != 0) 

table(dtf10$cluster_id)
cluster_dedup <- clusterfilt %>%
  dplyr::arrange(desc(Freq)) %>%
  distinct(Var2, .keep_all = T) %>%
  dplyr::arrange(Var1)

cluster_dedup
write_csv(cluster_dedup, "15aa_10groups_which_aa_per_cluster_distinct.csv")

# clustord <- clusterfilt[order(clusterfilt$Var1),]
# write_csv(clustord, "10groups_which_aa_per_cluster.csv")

# meds <- data.frame(cbind(1:15, pm$medoids, word(pm$medoids, 2, sep = "_")))
# write_csv(data.frame(meds), "output/cluster_key.csv")

pdf("output/15aa_10medoids_cluster_id_matrix.pdf", height = 20, width = 7)
ggplot(clusterdf, aes(Var1, Var2)) +
  geom_tile(aes(fill = Freq)) + 
  geom_text(aes(label = round(Freq, 1))) +
  scale_x_discrete(position = "top") +
  xlab("Cluster ID #") +
  ylab("AA monomer") +
  scale_fill_gradient(low = "white", high = "red") 
  #theme(axis.text.x = element_text(angle = 90))
        #axis.title.x = element_blank())
dev.off()
#write_csv(data.frame(table(dtf20$cluster_id, dtf20$substrate)),
#                     "output/table_of_aa_cluster_ids.csv")
# write_csv(dtf15, "output/amino_acid_15_clustering.csv")


# Hierarchical clustering
hier_clust <- hclust(mdist, method = "complete")


hier_clust$labels <- word(hier_clust$labels, 2, sep = "_")
# plot(hier_clust, cex = .7, labels = cdat$labels)
# hier_clust$labels <- cdat$labels
avg_dend_obj <- as.dendrogram(hier_clust)
avg_col_dend <- color_branches(hier_clust, k = 10)
avg_labl_dend <-  color_labels(avg_col_dend, k = 10)
avg_labl_dend  <- set(avg_labl_dend, "labels_cex", 0.5)
pdf("output/aa34_hierarchical_clustering_results.pdf", width = 80)
plot(avg_labl_dend)
dev.off()
col_labels <- get_leaves_branches_col(avg_labl_dend)

