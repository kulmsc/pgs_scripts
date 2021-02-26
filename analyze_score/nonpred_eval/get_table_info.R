


i <- 1
list_ukb_cc <- list()
list_gwas_cc <- list()

all_author <- unlist(lapply(strsplit(list.files("../tune_score/tune_results/", "res"), "_"), function(x) x[1]))
all_ss <- data.frame(matrix(0, nrow = length(all_author), ncol = 5))

for(author in all_author){

  ukb_pheno <- read.table(paste0("../construct_defs/pheno_defs/diag.", tolower(author) ,".txt.gz"), stringsAsFactors=F)

  pheno_eids <- read.table("../construct_defs/eid.csv", header = T)
  pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
  pheno_eids <- pheno_eids[-length(pheno_eids)]

  all_eid <- read.table("~/athena/doc_score/do_score/all_specs/for_eid.fam", stringsAsFactors=F)

  ukb_pheno <- ukb_pheno[pheno_eids %in% all_eid[,1],]

  all_ss[i, 2:3] <- as.numeric(table(rowSums(ukb_pheno[,1:4]) > 0))

  ss_cc <- read.table(paste0("~/athena/doc_score/raw_ss/", author, "/notes"), sep = "\t", stringsAsFactors = F)
  all_ss[i,4:5] <- as.numeric(ss_cc[4:5,])

  i <- i + 1
}

all_ss[,1] <- all_author
write.table(all_ss, "all_sample_size.txt", row.names = F, col.names = F, sep = " ", quote = F)
