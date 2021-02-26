library(survival)
library(data.table)
library(pROC)
library(epitools)

#author <- "liu-2"
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


all_auc <- rep(list(matrix(0, nrow = ncol(brit_scores), ncol = 6)), 8)
all_pval <- rep(list(matrix(0, nrow = ncol(brit_scores), ncol = 1)), 8)

df <- data.frame(phen = diag_phen, brit_covars)

census_attrib <- read.table("census_attrib.txt", stringsAsFactors=F, header=T)
colnames(census_attrib) <- c("eid", "median_age", "unemployed", "very_good_health", "pop_dens")
survey_attrib <- read.table("survey_attrib.txt", stringsAsFactors=F, header=T, sep=",")
colnames(survey_attrib) <- c("eid", "time_address", "number_in_house", "income", "age_edu")
common_eid <- intersect(df$eid[df$eid %in% survey_attrib$eid], df$eid[df$eid %in% census_attrib$eid])


brit_scores <- brit_scores[df$eid %in% common_eid,]
df <- df[df$eid %in% common_eid,]
survey_attrib <- survey_attrib[survey_attrib$eid %in% common_eid,]
census_attrib <- census_attrib[census_attrib$eid %in% common_eid,]
survey_attrib <- survey_attrib[order(survey_attrib$eid)[rank(df$eid)],]
census_attrib <- census_attrib[order(census_attrib$eid)[rank(df$eid)],]
df <- cbind(df, survey_attrib[,-1], census_attrib[,-1])

#QC #################################
df$time_address[df$time_address < 0 & !is.na(df$time_address)] <- NA
df$bin_address <- NA
df$bin_address[df$time_address %in% 1:19 & !is.na(df$time_address)] <- 0
df$bin_address[df$time_address %in% 20:100 & !is.na(df$time_address)] <- 1

df$income[(df$income == -1 | df$income == -3) & !is.na(df$income)] <- NA
df$bin_income <- NA
df$bin_income[df$income %in% c(1,2) & !is.na(df$income)] <- 0
df$bin_income[df$income %in% c(3,4,5) & !is.na(df$income)] <- 1

df$number_in_house[(df$number_in_house == -1 | df$number_in_house == -3) & !is.na(df$number_in_house)] <- NA
df$number_in_house[df$number_in_house > 25 & !is.na(df$number_in_house)] <- NA
df$bin_in_house <- NA
df$bin_in_house[df$number_in_house %in% c(1,2) & !is.na(df$number_in_house)] <- 0
df$bin_in_house[df$number_in_house %in% 3:25 & !is.na(df$number_in_house)] <- 1

df$age_edu[df$age_edu < 0 & !is.na(df$age_edu)] <- NA
df$bin_edu <- NA
df$bin_edu[df$age_edu %in% 1:19 & !is.na(df$age_edu)] <- 0
df$bin_edu[df$age_edu %in% 20:40 & !is.na(df$age_edu)] <- 1

df$bin_census_age <- NA
df$bin_census_age[df$median_age > median(df$median_age, na.rm = T) & !is.na(df$median_age)] <- 0
df$bin_census_age[df$median_age <= median(df$median_age, na.rm = T) & !is.na(df$median_age)] <- 1

df$bin_census_employ <- NA
df$bin_census_employ[df$unemployed > median(df$unemployed, na.rm = T) & !is.na(df$unemployed)] <- 0
df$bin_census_employ[df$unemployed <= median(df$unemployed, na.rm = T) & !is.na(df$unemployed)] <- 1

df$bin_census_health <- NA
df$bin_census_health[df$very_good_health > median(df$very_good_health, na.rm = T) & !is.na(df$very_good_health)] <- 0
df$bin_census_health[df$very_good_health <= median(df$very_good_health, na.rm = T) & !is.na(df$very_good_health)] <- 1

df$bin_pop_den <- NA
df$bin_pop_den[df$pop_dens > median(df$pop_dens, na.rm = T) & !is.na(df$pop_dens)] <- 0
df$bin_pop_den[df$pop_dens <= median(df$pop_dens, na.rm = T) & !is.na(df$pop_dens)] <- 1
#####################################



vars_examine <- grep("bin", colnames(df), value=T)

for(j in 1:ncol(brit_scores)){
  df$score <- brit_scores[,j]
  train_df <- df[df$eid %in% brit_train[,1],]
  test_df <- df[df$eid %in% brit_test[,1],]

  for(var in vars_examine){
    for(bin_val in 0:1){

    sub_train_df <- train_df[train_df[var] == bin_val,]
    sub_test_df <- test_df[test_df[var]== bin_val,]

    if("sex" %in% colnames(sub_train_df)){
      mod <- glm(phen ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = sub_train_df, family = "binomial")
    } else {
      mod <- glm(phen ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = sub_train_df, family = "binomial")
    }
    
    test_roc <- roc(sub_test_df$phen ~ predict(mod, sub_test_df))

    if(bin_val == 0){
      save_roc <- test_roc
    } else {
      all_pval[[which(vars_examine == var)]][j,1] <- roc.test(save_roc, test_roc)$p.value
    }

    all_auc[[which(vars_examine == var)]][j,(bin_val*3+1):(bin_val*3+3)] <- as.numeric(ci.auc(test_roc))

    }
  }
}


###########################

so_far_done <- readRDS(paste0("per_score_results/", tolower(author), ".res.ethnictest.RDS"))
so_far_done[["many_auc"]] <- all_auc
so_far_done[["many_pval"]] <- all_pval
saveRDS(so_far_done, paste0("per_score_results/", tolower(author), ".res.many.RDS"))
