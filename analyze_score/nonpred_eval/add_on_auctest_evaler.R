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
escores <- apply(escores, 2, function(x) (x-mean(x)) / (max(abs((x-mean(x)))) * 0.01) )
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
########### SEX DIFFERENCES ####################
##############################################

diag_phen <- apply(brit_pheno[,1:5], 1, function(x) any(x==1))*1

if(sex_group == "A"){
  all_sex_auc <- matrix(0, nrow = ncol(brit_scores), ncol = 6)

  df <- data.frame(phen = diag_phen, brit_covars)
  sex_auc_p <- rep(NA, ncol(brit_scores))
  for(j in 1:ncol(brit_scores)){
    df$score <- brit_scores[,j]
    train_df <- df[df$eid %in% brit_train[,1],]
    test_df <- df[df$eid %in% brit_test[,1],]
    roc_holder <- list()
    for(sex in c(0,1)){
      sub_train_df <- train_df[train_df$sex == sex,]
      sub_test_df <- test_df[test_df$sex == sex,]
      sub_train_df <- sub_train_df[,-which(colnames(sub_train_df) == "sex")]
      sub_test_df <- sub_test_df[,-which(colnames(sub_test_df) == "sex")]
      mod <- glm(phen ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = sub_train_df, family = "binomial")
      test_roc <- roc(sub_test_df$phen ~ predict(mod, sub_test_df))
      all_sex_auc[j,(sex*3+1):(sex*3+3)] <- as.numeric(ci.auc(test_roc))
      roc_holder[[sex+1]] <- test_roc
    }
    sex_auc_p[j] <- roc.test(roc_holder[[1]], roc_holder[[2]], paired = FALSE)$p.value
  }
} else {
  all_sex_auc <- NULL
  sex_auc_p <- NULL
}

return_sex_split <- list("all_sex_auc" = all_sex_auc)


##############################################
########### AGE DIFFERENCE ####################
##############################################

split_age <- function(indf){
  quota <- sum(df$phen)/2
  age_cut <- min(indf$age)
  good_val <- sum(indf$phen[indf$age <= age_cut])
  while(good_val < quota){
    age_cut <- age_cut + 0.1
    good_val <- sum(indf$phen[indf$age <= age_cut])
  }
  indf$young_old <- 0
  indf$young_old[indf$age > age_cut] <- 1
  return(indf)
}

diag_phen <- apply(brit_pheno[,1:5], 1, function(x) any(x==1))*1

all_age_auc <- matrix(0, nrow = ncol(brit_scores), ncol = 6)

df <- data.frame(phen = diag_phen, brit_covars)
df <- split_age(df)
age_auc_p <- rep(NA, ncol(brit_scores))
for(j in 1:ncol(brit_scores)){
  df$score <- brit_scores[,j]
  train_df <- df[df$eid %in% brit_train[,1],]
  test_df <- df[df$eid %in% brit_test[,1],]
  roc_holder <- list()
  for(young in c(0,1)){
    sub_train_df <- train_df[train_df$young_old == young,]
    sub_test_df <- test_df[test_df$young_old == young,]
    sub_train_df <- sub_train_df[,-which(colnames(sub_train_df) == "young_old")]
    sub_test_df <- sub_test_df[,-which(colnames(sub_test_df) == "young_old")]
    if("sex" %in% colnames(sub_train_df)){
      mod <- glm(phen ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = sub_train_df, family = "binomial")
    } else {
      mod <- glm(phen ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = sub_train_df, family = "binomial")
    }
    test_roc <- roc(sub_test_df$phen ~ predict(mod, sub_test_df))
    all_age_auc[j,(young*3+1):(young*3+3)] <- as.numeric(ci.auc(test_roc))
    roc_holder[[young+1]] <- test_roc
  }
  age_auc_p[j] <- roc.test(roc_holder[[1]], roc_holder[[2]], paired = FALSE)$p.value
}




###########################

so_far_done <- readRDS(paste0("per_score_results/", tolower(author), ".res.wage.RDS"))
so_far_done[["auc_p"]] <- list("sex" = sex_auc_p, "age" = age_auc_p)
saveRDS(so_far_done, paste0("per_score_results/", tolower(author), ".res.auctest.RDS"))
