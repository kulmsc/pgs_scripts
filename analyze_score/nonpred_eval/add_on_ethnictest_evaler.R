library(survival)
library(data.table)
library(pROC)
library(epitools)

#author <- "Malik"
author <- commandArgs(trailingOnly=TRUE)

#read in scores
all_scores <-  readRDS(paste0("../../do_score/final_scores/all_score.", tolower(author), ".RDS"))
score_eids <- read.table("~/athena/doc_score/do_score/all_specs/for_eid.fam", stringsAsFactors=F)

other_scores <- readRDS(paste0("~/athena/doc_score/do_score/ethnic_scoring/final_scores/all_score.", tolower(author), ".RDS"))
ethnic_eids <- read.table("~/athena/doc_score/qc/ethnic_eids", stringsAsFactors=F)

#normalize the scores
scores <- all_scores[,grepl(tolower(author), colnames(all_scores))]
scores <- scores[,apply(scores, 2, function(x) length(unique(x)) > 3)]
scores <- apply(scores, 2, function(x) (x-mean(x)) / (max(abs((x-mean(x)))) * 0.01) )
rm(all_scores)

escores <- other_scores[,grepl(tolower(author), colnames(other_scores))]
escores <- escores[,apply(escores, 2, function(x) length(unique(x)) > 3)]
#escores <- apply(escores, 2, function(x) (x-mean(x)) / (max(abs((x-mean(x)))) * 0.01) )
rm(other_scores)

#read in phenotype and covars
pheno <- read.table(paste0("../construct_defs/pheno_defs/diag.", tolower(author), ".txt.gz"), stringsAsFactors=F)
dates <- read.table(paste0("../construct_defs/pheno_defs/time.", tolower(author), ".txt.gz"), stringsAsFactors=F)
covars <- readRDS("../get_covars/base_covars.RDS")

#sort the phenotypes
pheno_eids <- read.table("../construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]

#split covars
brit_covars <- covars[(covars[,1] %in% score_eids[,1]) & (covars[,1] %in% pheno_eids),]
brit_covars <- brit_covars[order(brit_covars[,1]),]

ethnic_covars <- covars[(covars[,1] %in% ethnic_eids[,1]) & (covars[,1] %in% pheno_eids),]
ethnic_covars <- ethnic_covars[order(ethnic_covars[,1]),]

#split pheno, dates
brit_pheno <- pheno[(pheno_eids %in% score_eids[,1]) & (pheno_eids %in% covars[,1]),]
brit_dates <- dates[(pheno_eids %in% score_eids[,1]) & (pheno_eids %in% covars[,1]),]
brit_pheno_eids <- pheno_eids[(pheno_eids %in% score_eids[,1]) & (pheno_eids %in% covars[,1])]

ethnic_pheno <- pheno[(pheno_eids %in% ethnic_eids[,1]) & (pheno_eids %in% covars[,1]),]
ethnic_dates <- dates[(pheno_eids %in% ethnic_eids[,1]) & (pheno_eids %in% covars[,1]),]
ethnic_pheno_eids <- pheno_eids[(pheno_eids %in% ethnic_eids[,1]) & (pheno_eids %in% covars[,1])]

#split scores
brit_scores <- scores[(score_eids[,1] %in% pheno_eids) & (score_eids[,1] %in% covars[,1]),]
brit_score_eids <- score_eids[(score_eids[,1] %in% pheno_eids) & (score_eids[,1] %in% covars[,1]),]

brit_scores <- brit_scores[order(brit_score_eids[,1]),]
brit_eid <- brit_score_eids[order(brit_score_eids[, 1]), 1]

ethnic_scores <- escores[(ethnic_eids[,1] %in% pheno_eids) & (ethnic_eids[,1] %in% covars[,1]),]
ethnic_score_eids <- ethnic_eids[(ethnic_eids[,1] %in% pheno_eids) & (ethnic_eids[,1] %in% covars[,1]),,drop=F]

ethnic_scores <- ethnic_scores[order(ethnic_score_eids[,1]),]
ethnic_eid <- ethnic_score_eids[order(ethnic_score_eids[, 1]), 1]

rm(scores)
rm(escores)

#get the british eids
brit_train <- read.table("../../qc/cv_files/train_eid.0.6.txt", stringsAsFactors=F)
brit_test <- read.table("../../qc/cv_files/test_eid.0.4.txt", stringsAsFactors=F)

#remove possible sex group
author_defs <- read.table("../descript_defs/author_defs", stringsAsFactors=F, header=T)
sex_group <- author_defs[author_defs[,1] == author, 3]

if(sex_group == "M"){
  all_sex <- read.table("all_sex", stringsAsFactors=F, sep=",", header=T)
  all_sex <- all_sex[all_sex[,1] %in% ethnic_eid,]
  all_sex <- all_sex[order(all_sex[,1])[rank(ethnic_eid)],]

  brit_scores <- brit_scores[brit_covars$sex == 1,]
  brit_eid <- brit_eid[brit_covars$sex == 1]
  brit_pheno <- brit_pheno[brit_covars$sex == 1,]
  brit_dates <- brit_dates[brit_covars$sex == 1,]
  brit_covars <- brit_covars[brit_covars$sex == 1,]

  ethnic_scores <- ethnic_scores[all_sex[,2] == 1,]
  ethnic_eid <- ethnic_eid[all_sex[,2] == 1]
  ethnic_pheno <- ethnic_pheno[all_sex[,2] == 1,]
  ethnic_dates <- ethnic_dates[all_sex[,2] == 1,]
  ethnic_covars <- ethnic_covars[all_sex[,2] == 1,]

}else if(sex_group == "F"){
  all_sex <- read.table("all_sex", stringsAsFactors=F, sep=",", header=T)
  all_sex <- all_sex[all_sex[,1] %in% ethnic_eid,]
  all_sex <- all_sex[order(all_sex[,1])[rank(ethnic_eid)],]

  brit_scores <- brit_scores[brit_covars$sex == 0,]
  brit_eid <- brit_eid[brit_covars$sex == 0]
  brit_pheno <- brit_pheno[brit_covars$sex == 0,]
  brit_dates <- brit_dates[brit_covars$sex == 0,]
  brit_covars <- brit_covars[brit_covars$sex == 0,]

  ethnic_scores <- ethnic_scores[all_sex[,2] == 0,]
  ethnic_eid <- ethnic_eid[all_sex[,2] == 0]
  ethnic_pheno <- ethnic_pheno[all_sex[,2] == 0,]
  ethnic_dates <- ethnic_dates[all_sex[,2] == 0,]
  ethnic_covars <- ethnic_covars[all_sex[,2] == 0,]
}



##############################################
########### RACIAL GROUPS ####################
##############################################


asian_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/asian_fam", stringsAsFactors=F)
african_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/african_fam", stringsAsFactors=F)
euro_fam <- read.table("~/athena/ukbiobank/custom_qc/fam_files/euro_fam", stringsAsFactors=F)

ethnic_scores <- as.data.frame(ethnic_scores)
ethnic_scores$ethnic <- "none"
ethnic_scores$ethnic[ethnic_eid %in% asian_fam[,1]] <- "asian"
ethnic_scores$ethnic[ethnic_eid %in% african_fam[,1]] <- "african"
ethnic_scores$ethnic[ethnic_eid %in% euro_fam[,1]] <- "euro"

normalize <- function(x){
 (x - min(x))/((max(x) - min(x)) * 0.01)
}

for(i in 1:(ncol(ethnic_scores)-1)){
  ethnic_scores[ethnic_scores$ethnic == "asian", i] <- normalize(ethnic_scores[ethnic_scores$ethnic == "asian", i])
  ethnic_scores[ethnic_scores$ethnic == "african", i] <- normalize(ethnic_scores[ethnic_scores$ethnic == "african", i])
  ethnic_scores[ethnic_scores$ethnic == "euro", i] <- normalize(ethnic_scores[ethnic_scores$ethnic == "euro", i])
}

for(i in 1:ncol(brit_scores)){
  brit_scores[,i] <- normalize(brit_scores[,i])
}


brit_stats <- matrix(0, nrow = ncol(ethnic_scores)-1, ncol = 13)
euro_stats <- matrix(0, nrow = ncol(ethnic_scores)-1, ncol = 13)
african_stats <- matrix(0, nrow = ncol(ethnic_scores)-1, ncol = 13)
asian_stats <- matrix(0, nrow = ncol(ethnic_scores)-1, ncol = 13)

for(i in 1:(ncol(ethnic_scores)-1)){
  brit_ind <- which(colnames(brit_scores) == colnames(ethnic_scores)[1])
  brit_stats[i,] <- c(mean(brit_scores[,brit_ind]), sd(brit_scores[,brit_ind]), quantile(brit_scores[,brit_ind], 0:10/10))
  sub_score <- ethnic_scores[ethnic_scores$ethnic == "euro",i]
  euro_stats[i,] <- c(mean(sub_score), sd(sub_score), quantile(sub_score, 0:10/10))
  sub_score <- ethnic_scores[ethnic_scores$ethnic == "african",i]
  african_stats[i,] <- c(mean(sub_score), sd(sub_score), quantile(sub_score, 0:10/10))
  sub_score <- ethnic_scores[ethnic_scores$ethnic == "asian",i]
  asian_stats[i,] <- c(mean(sub_score), sd(sub_score), quantile(sub_score, 0:10/10))
}


tvals <- data.frame(matrix(NA, nrow = 6, ncol = ncol(ethnic_scores)-1))
for(i in 1:(ncol(ethnic_scores)-1)){
  brit_ind <- which(colnames(brit_scores) == colnames(ethnic_scores)[i])
  euro <- ethnic_scores[ethnic_scores$ethnic == "euro",i]
  african <- ethnic_scores[ethnic_scores$ethnic == "african",i]
  asian <- ethnic_scores[ethnic_scores$ethnic == "asian",i]
  brit <- brit_scores[,brit_ind]

  tvals[1,i] <- t.test(brit, euro)$stat
  tvals[2,i] <- t.test(brit, african)$stat
  tvals[3,i] <- t.test(brit, asian)$stat
  tvals[4,i] <- t.test(euro, asian)$stat
  tvals[5,i] <- t.test(euro, african)$stat
  tvals[6,i] <- t.test(african, asian)$stat
}

pvals <- data.frame(matrix(NA, nrow = 6, ncol = ncol(ethnic_scores)-1))
for(i in 1:(ncol(ethnic_scores)-1)){
  brit_ind <- which(colnames(brit_scores) == colnames(ethnic_scores)[i])
  euro <- ethnic_scores[ethnic_scores$ethnic == "euro",i]
  african <- ethnic_scores[ethnic_scores$ethnic == "african",i]
  asian <- ethnic_scores[ethnic_scores$ethnic == "asian",i]
  brit <- brit_scores[,brit_ind]

  pvals[1,i] <- t.test(brit, euro)$p.value
  pvals[2,i] <- t.test(brit, african)$p.value
  pvals[3,i] <- t.test(brit, asian)$p.value
  pvals[4,i] <- t.test(euro, asian)$p.value
  pvals[5,i] <- t.test(euro, african)$p.value
  pvals[6,i] <- t.test(african, asian)$p.value
}



tvals$comp1 <- c("brit", "brit", "brit", "euro", "euro", "african")
tvals$comp2 <- c("euro", "african", "asian", "asian", "african", "asian")

return_ethnic_groups <- list("names_to_keep" = colnames(ethnic_scores), "brit_stats" = brit_stats, "euro_stats" = euro_stats, "african_stats" = african_stats, "asian_stats" = asian_stats, "stat_test" = tvals, "pval_test" = pvals)





###########################

so_far_done <- readRDS(paste0("per_score_results/", tolower(author), ".res.auctest.RDS"))
so_far_done[["new_ethnic"]] <- return_ethnic_groups
saveRDS(so_far_done, paste0("per_score_results/", tolower(author), ".res.ethnictest.RDS"))
