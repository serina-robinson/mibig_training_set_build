## Install packages
pacman::p_load('ape', 'dendextend', 'rgl', 'cluster', 'ggrepel', 'plotly', 'tidyverse', 'data.table', 'readxl', 'stringr', 'Biostrings', 'DECIPHER', 'phangorn', 'ggplot2', 'seqinr', 'bgafun')

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read in the dataset
dat <- read_delim("data/sp2_34extract_named_features.txt", delim = "\t", col_names = F)
head(dat)

# Numerical data
mat <- as.matrix(dat[,3:ncol(dat)])
rownames(mat) <- dat$X2


# Principal components
dat.pca <- prcomp(mat, center = TRUE, scale = TRUE)
summary(dat.pca)
# 
# # library(devtools)
# # install_github("vqv/ggbiplot")
# library(ggbiplot)
ggbiplot(dat.pca, labels = NULL, loadings = F, var.axes = F, loadings.label = FALSE)
ggscreeplot(dat.pca)
# plot(dat.pca$x[,47], dat.pca$x[,48])
# Physicochemical dist
 #mdist <- dist(mat)

x<-dat.pca$x[,1]
y<-dat.pca$x[,2]
z<-dat.pca$x[,3]

# Calculate percent explained by each principal component
eigs <- dat.pca$sdev^2

############# 
# Do sequences cluster by their substrate?
dtf <- data.frame(dat$X1, dat$X2, stringsAsFactors = F)
head(dtf)
colnames(dtf) <- c("nms", "subst")

nrow(dtf)
dat <- data.frame(cbind(x, y, dtf$subst), stringsAsFactors = F)
# rownames(dat) <- rd$X1
rownames(dat) <- dtf$nms
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
cdat$labels

pl <- list()
for(i in 1:5) {
  pl[[i]] <- list()
  for(j in 1:5) {
    cdat$x <- dat.pca$x[,i]
    cdat$y <- dat.pca$x[,j]
    pc1 <- eigs[i]/sum(eigs)
    pc2 <- eigs[j]/sum(eigs)
    xlab<-paste0("PC ", i, " (", round( pc1 * 100, 2), "% explained var.)")
    ylab<-paste0("PC ", j, " (", round( pc2 * 100, 2), "% explained var.)")
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


pdf("data/15properties_PCA_condensed.pdf")
pl[[1]][[2]]
dev.off()

pdf("data/15properties_PCA_aa_plot.pdf", width = 50, height = 50)
plot_grid(pl[[1]][[1]], pl[[1]][[2]], pl[[1]][[3]], pl[[1]][[4]], pl[[1]][[5]],
          pl[[2]][[1]], pl[[2]][[2]], pl[[2]][[3]], pl[[2]][[4]], pl[[2]][[5]],
          pl[[3]][[1]], pl[[3]][[2]], pl[[3]][[3]], pl[[3]][[4]], pl[[3]][[5]],
          pl[[4]][[1]], pl[[4]][[2]], pl[[4]][[3]], pl[[4]][[4]], pl[[4]][[5]],
          pl[[5]][[1]], pl[[5]][[2]], pl[[5]][[3]], pl[[5]][[4]], pl[[5]][[5]], 
          nrow = 5, ncol = 5)
dev.off()

# PAM medoids (Partitioning around medoids)

pm <- pam(mat, k = 10, metric = "euclidean", stand = FALSE)

dtf10 <- data.frame(rownames(cdat), names(pm$clustering),  pm$clustering)
colnames(dtf10) <- c("full_id", "substrate", "cluster_id")
head(dtf10)
# write_csv(dtf10, "output/amino_acid_10_clustering.csv")

pm <- pam(mat, k = 6, metric = "euclidean", stand = FALSE)
dtf6 <- data.frame(rownames(cdat), names(pm$clustering),  pm$clustering)
colnames(dtf6) <- c("full_id", "substrate", "cluster_id")
table(dtf6$cluster_id)

# write_csv(dtf15, "output/amino_acid_15_clustering.csv")

# pm <- pam(mat, k = 16, metric = "euclidean", stand = FALSE)
# dtf16 <- data.frame(rownames(cdat), names(pm$clustering),  pm$clustering)
# colnames(dtf16) <- c("full_id", "substrate", "cluster_id")
# table(dtf16$cluster_id)
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
table(dtf6$cluster_id)
dtftabs[grep("iso", dtftabs$substrate),]
dtftabs <- dtf6 %>%
  group_by(cluster_id, substrate) %>%
  add_count(substrate) 
table(dtftabs)  

clusterdf <- data.frame(table(dtf6$cluster_id, dtf12$substrate))
clusterfilt <- clusterdf %>%
  dplyr::filter(Freq != 0)

colnames(clusterfilt) <- c("cluster_id", "aa_monomer", "frequency")
# clustersort <- clusterfilt %>%
#   spread(., cluster_id, value = frequency)

# Col sums
col_sums <- colSums(clustersort[,2:16])
col_sums
col_maxs <- apply(clustersort[,2:16], 2, max)
col_sums - col_maxs
col_maxs
round(col_maxs/(col_sums-col_maxs), 2)
# Row sums
row_sums <- rowSums(clustersort[,2:16])
row_sums
row_maxs <- apply(clustersort[,2:16], 1, max)
row_maxs
row_sums - row_maxs

colnames(clusterfilt) <- c("cluster_id", "aa_monomer", "frequency")
cluster_dedup <- clusterfilt %>%
  group_by(cluster_id) %>%
  do( data.frame(with(data=., .[order(frequency, decreasing = T),] )) )
  # spread(., cluster_id, value = frequency)
tail(cluster_dedup)

# dplyr::arrange(desc(Freq)) %>%
  # distinct(Var2, .keep_all = T) %>%
  # dplyr::arrange(Var1)
head(cluster_dedup)
write_csv(cluster_dedup, "15aa_15groups_all_included_aa_per_cluster_distinct.csv")

# clustord <- clusterfilt[order(clusterfilt$Var1),]
# write_csv(clustord, "10groups_which_aa_per_cluster.csv")

# meds <- data.frame(cbind(1:15, pm$medoids, word(pm$medoids, 2, sep = "_")))
# write_csv(data.frame(meds), "output/cluster_key.csv")

pdf("output/15aa_6medoids_cluster_id_matrix.pdf", height = 20, width = 7)
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
# with(cdat, plot3d(x = x, y = y, z = z,
#                   col = colrs, type = 's', size = 1))


