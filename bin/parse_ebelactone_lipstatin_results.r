# Install packages
pacman::p_load('data.table', 'Biostrings', 'DECIPHER', 'rentrez', 'genbankr', 'RCurl', 'tidyverse', 'XML', 'RDSTK', 'rvest', 'ape')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("genbankr", version = "3.8")
setwd("~/Documents/Wageningen_UR/github/mibig_training_set_build_test/")

ebe <- fread("~/Downloads/ebelactone_asdb_search_results.csv", data.table =F, sep = "\t", quote = "") %>%
  janitor::clean_names() %>%
  mutate(download_url = gsub('\\"', "", download_url))

lip <- fread("~/Downloads/lipstatin_asdb_search_results.csv", data.table =F, sep = "\t", quote = "") %>%
  janitor::clean_names() %>%
  mutate(download_url = gsub('\\"', "", download_url))

bnd <- ebe %>%
  bind_rows(lip)
bnd$ncbi_accession
genbank_acc <- word(word(bnd$download_url, 1, sep = "cluster"), sep = "genbank/", 2)

a <- bnd$cluster_number
bnd$download_url
genbank_fils <- paste0("https://antismash-db.secondarymetabolites.org/output/", genbank_acc, bnd$ncbi_accession, ".", sprintf("cluster%03d", a), ".gbk")
genbank_fils[1] <- 'https://antismash-db.secondarymetabolites.org/output/GCF_000010725/NC_013854.1.cluster003.gbk'
shell_commands <- paste0("curl -o ", bnd$ncbi_accession, " ", bnd$download_url)
shell_commands
write(shell_commands, "output/genbank_fil_paths.sh")
shell_commands

# Read in the genbank HTML
gbhtml <- list.files("output/antismashdb2_test_fil_paths", full.names = T)

scraping_wiki <- gbhtml %>%
  map(read_html) %>%
  purrr::map(~ html_nodes(.x, "a")) %>%
  purrr::map(~ html_attr(.x, "href")) %>%
  unlist()

gbk_to_pull <- paste0('https://antismash-db.secondarymetabolites.org/', scraping_wiki)
gbk_to_pull
gbk_sh <- write(paste0("curl -O ", gbk_to_pull), "output/genbank_fil_path_final.sh")


#   map(~html_attr("href", x))

# Pull the genbank files for all 
# gnbnk <- paste0(ebe$download_url[1], ebe$ncbi_accession[1])
# tmp <- getURL('https://antismash-db.secondarymetabolites.org/output/GCF_000010725/NC_013854.1.cluster003.gbk', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)


gb <- list.files("output/antismashdb2_lip_ebe_pull", full.names =  T)

# trans <- gb_sample$FEATURES[grep("translation", gb_sample$FEATURES)]
gb2 <- gb[c(-146, -182)]
parseGenBank(gb[200])

dtf <- list()
# for(i in 1:length(gb)) {
for(i in 1:length(gb2)) {
  prs <- parseGenBank(gb2[i])
  feats <- prs$FEATURES[grep("translation", prs$FEATURES)]
  dtf[[i]] <- list()
  for(j in 1:length(feats)) {
    try ({
    tobind <- data.frame(feats[j]) %>%
      dplyr::select(#contains("translation"),
                    contains("protein_id"),
                    contains("product"))
    dtf[[i]][[j]] <- tobind
    })
  }
}

# Unlist everything
dtf2 <- dtf
dtf3 <- unlist(dtf2)
# dtf3
# dtf3
# head(dtf3)
# inds_topull <- grep("peptide synthase|peptide synthetase|CoA ligase", dtf3) + 1
# 
# amp <- names(dtf3)[grep("peptide synthase|peptide synthetase|CoA ligase", dtf3)]
# sum(table(dtf3[grep("peptide synthase|peptide synthetase|CoA ligase", dtf3)]))
# length(dtf3[grep("peptide synthase|peptide synthetase|CoA ligase", dtf3)])
# 
# 
# nums <- paste0("X", parse_number(amp), ".protein_id")
# nums
# topull <- as.character(dtf3[inds_topull])
# 
# write(topull, "output/dtf_accs_for_entrez.txt")
# 
# # dtfmat <- data.frame(matrix(unlist(dtf), ncol=2),stringsAsFactors=FALSE)


dtfunl <- data.frame(unlist(dtf), stringsAsFactors = F) %>%
  dplyr::filter(grepl("\\.1", unlist.dtf.))


dtf_accs <- unique(as.character(dtfunl[,1]))
length(dtf_accs)

# write(dtf_accs, "output/dtf_accs_for_entrez.txt")
# 
# test_acc <- entrez_fetch(topull,
#               db = "protein", rettype = "fasta", ap_key="826357f5ff17c7ec62e583909071e94f9d08")
# 

# Find intersection between BLAST hits for LstC and orf1 and the results
lstc_orf1 <- read_csv("~/Downloads/973Y688E014-Alignment-HitTable.csv", col_names = F)
lst_orf <- dtf_accs[dtf_accs %in% lstc_orf1$X2]

efetchr <- entrez_fetch(lst_orf,
db = "protein", rettype = "fasta", ap_key="826357f5ff17c7ec62e583909071e94f9d08")
write(efetchr, "output/lstc_orf_antismashdb_hits.fasta")
faaa <- readAAStringSet("output/lstc_orf_antismashdb_hits.fasta")
names(faaa)

# Ebelactone
ebe <- ebe %>%
  dplyr::filter(most_similar_known_cluster == "Ebelactone") 

ebe_sort <- ebe[order(ebe$similarity_in_percent, decreasing = T),]
head(ebe_sort)

lip <- lip %>%
  janitor::clean_names() %>%
  dplyr::filter(most_similar_known_cluster == "Lipstatin") 

lip_sort <- lip[order(lip$similarity_in_percent, decreasing = T),] 

to_write <- ebe_sort %>%
  bind_rows(lip_sort) 

to_write$number_genus <- gsub('\\"', "", to_write$number_genus)
write_csv(to_write, "output/top_ebelactone_lipstatin_hits.csv")

