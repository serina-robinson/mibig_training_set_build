color_pcoa_points <- function(namvec, dat, pal) {
  # namvec must be a partial match in the rownames of the pcoa data

  colrs <- rep("gray", nrow(dat))
  alph  <- rep(0.75, nrow(dat))
  sz <- rep(1, nrow(dat))
  nms <- rep("unknown", nrow(dat))
  
  for(i in 1:length(namvec)) {
    inds <- grep(namvec[i], rownames(dat))
    colrs[inds] <- pal[i]
    alph[inds] <- 1
    # sz[inds] <- 2
    nms[inds] <- namvec[i]
  }
  
  dat$colrs <- colrs
  dat$alph <- alph
  dat$sz <- sz
  dat$nms <- nms
  return(dat)
}
