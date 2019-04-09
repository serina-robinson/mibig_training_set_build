# Read in the reference sequence
query <- readAAStringSet("../data/test1.faa") 
ref <- readAAStringSet("../data/A_domains_muscle.fasta")

# Align the query and the refereuence
# alned <- AlignSeqs(c(query, ref))
alned <- AAStringSet(muscle(c(query, ref), in1 = "../data/A_domains_muscle.fasta", in2 = "../data/test1.faa", profile = T))
BrowseSeqs(alned)
names(alned)
query_aln <- alned[1]
ref_aln <- alned["P0C062_A1"]

# Read in the 34 indices
aa34 <- fread("../data/A34positions.txt", data.table = F, header = F)
aa34_inds <- as.numeric(aa34[1,])
aa34_inds_adj <- aa34_inds - 65

# Exract the 34 amino acid positions
poslist <- list()
position = 1

for(i in 1:width(ref_aln)) {
  if (substr(ref_aln, i, i) != "-") {
    if (position %in% aa34_inds_adj) {
      poslist[[i]] <- i
    }
    position = position + 1
  }
}

# Get the new indices
new_34inds <- unlist(poslist)

# Get 34 aa code
query_pos <- as.character(unlist(lapply(1:length(new_34inds), function(x) {
  substr(query_aln, new_34inds[x], new_34inds[x]) })))
paste0(query_pos, collapse = "")

feats <- convert_seq_15aap(query_pos)

# Center and scale
feats_scaled <- scale(feats, center = T)
names(feats_scaled) <- colnames(feats)
return(feats_scaled)