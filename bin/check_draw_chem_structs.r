# Install packages
pacman::p_load("webchem", "readxl", "tidyverse", "ChemmineR", "ChemmineOB", "data.table", "plot3D", "plot3Drgl",
               "ggrepel")

# Set working directory
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

# Read SDFset
sdfset <- read.SDFset("data/90_cids_sdf.sdf")

# Plot the structure
sdfset2 <- regenerateCoords(sdfset[1:5])
pdf("output/chemical_structures_plotted.pdf", width = 50, height = 50)
plot(sdfset[1:90], print = F)
dev.off()
# openBabelPlot(sdfset[1],regenCoords=TRUE)
# BiocManager::install("ChemmineOB")
