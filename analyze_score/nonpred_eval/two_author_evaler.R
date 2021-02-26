library(survival)
library(data.table)
library(pROC)
library(epitools)

author <- "Michailidou"

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
########### SCORE SIZES ####################
##############################################

splitter <- function(values, N){
  inds = c(0, sapply(1:N, function(i) which.min(abs(cumsum(as.numeric(values)) - sum(as.numeric(values))/N*i))))
  dif = diff(inds)
  if(any(dif == 0)){
    dif[1] <- dif[1] - sum(dif == 0)
    dif[dif == 0] <- 1
  }

  re = rep(1:length(dif), times = dif)
  return(split(values, re))
}


all_equal_splits_len <- matrix(0, nrow = ncol(brit_scores), ncol = 4)
all_equal_splits_mean <- matrix(0, nrow = ncol(brit_scores), ncol = 4)
all_quantiles <- matrix(0, nrow = ncol(brit_scores), ncol = 10)
all_len <- rep(0, ncol(brit_scores))

for(i in 1:ncol(brit_scores)){
  score_name <- colnames(brit_scores)[i]
  split_name <- strsplit(score_name, ".", fixed = T)[[1]]
  score_files <- list.files(paste0("~/athena/doc_score/mod_sets/", author, "/"), glob2rx(paste0(split_name[1], "*", split_name[3], ".", split_name[2], ".ss")))
  score_list <- list()
  for(j in 1:length(score_files)){
    #score_list[[j]] <- read.table(paste0("~/athena/doc_score/mod_sets/", author, "/", score_files[j]), stringsAsFactors = F, header = F)
    score_list[[j]] <- as.data.frame(fread(paste0("~/athena/doc_score/mod_sets/", author, "/", score_files[j])))
    colnames(score_list[[j]]) <- c("CHR", "BP", "RSID", "A1",  "A2", "SE", "BETA", "P", "ESS")
  }

  ss <- do.call("rbind", score_list)
  if(ss[1,7] == "BETA"){
    ss <- ss[ss[,7] != "BETA",]
    ss[,7] <- as.numeric(ss[,7])
  }
  ss[,7] <- as.numeric(ss[,7])
  ss <- ss[!is.na(ss[,7]) & ss[,7] != 0,]
  ss$eff <- abs(ss[,7])
  ss <- ss[order(ss$eff),]

  all_equal_splits_len[i,] <- unlist(lapply(splitter(ss$eff, 4), length))
  all_equal_splits_mean[i,] <- unlist(lapply(splitter(ss$eff, 4), mean))
  all_quantiles[i,] <- quantile(ss$eff, 1:10/10)
  all_len[i] <- length(ss$eff)
}

return_score_sizes <- list("all_equal_splits_len" = all_equal_splits_len, "all_equal_splits_mean" = all_equal_splits_mean, "all_quantiles" = all_quantiles, "all_len" = all_len)



##############################################
########### PHENOTYPE DEFINITIONS ####################
##############################################

icd_phen <- apply(brit_pheno[,3:4], 1, function(x) any(x==1))*1
selfrep_phen <- apply(brit_pheno[,1:2], 1, function(x) any(x==1))*1
diag_phen <- apply(brit_pheno[,1:5], 1, function(x) any(x==1))*1
any_phen <- apply(brit_pheno, 1, function(x) any(x==1))*1
double_phen <- apply(brit_pheno, 1, function(x) sum(x) > 1)*1
all_phens <- list(icd_phen, selfrep_phen, diag_phen, any_phen, double_phen)

all_phen_auc <- matrix(0, nrow = ncol(brit_scores), ncol = 15)
phen_type_starts <- seq(1, 15, 3)

for(i in 1:5){
  df <- data.frame(phen = all_phens[[i]], brit_covars)
  for(j in 1:ncol(brit_scores)){
    df$score <- brit_scores[,j]
    train_df <- df[df$eid %in% brit_train[,1],]
    test_df <- df[df$eid %in% brit_test[,1],]
    mod <- glm(phen ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = train_df, family = "binomial")
    test_roc <- roc(test_df$phen ~ predict(mod, test_df))
    all_phen_auc[j,(phen_type_starts[i]):(phen_type_starts[i]+2)] <- as.numeric(ci.auc(test_roc))
  }
}

return_pheno_defs <- list("all_phen_auc" = all_phen_auc, "total_phens" = unlist(lapply(all_phens, sum)))


##############################################
########### SEX DIFFERENCES ####################
##############################################

if(sex_group == "A"){
  all_sex_auc <- matrix(0, nrow = ncol(brit_scores), ncol = 6)

  df <- data.frame(phen = diag_phen, brit_covars)
  for(j in 1:ncol(brit_scores)){
    df$score <- brit_scores[,j]
    train_df <- df[df$eid %in% brit_train[,1],]
    test_df <- df[df$eid %in% brit_test[,1],]
    for(sex in c(0,1)){
      sub_train_df <- train_df[train_df$sex == sex,]
      sub_test_df <- test_df[test_df$sex == sex,]
      sub_train_df <- sub_train_df[,-which(colnames(sub_train_df) == "sex")]
      sub_test_df <- sub_test_df[,-which(colnames(sub_test_df) == "sex")]
      mod <- glm(phen ~ age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + score, data = sub_train_df, family = "binomial")
      test_roc <- roc(sub_test_df$phen ~ predict(mod, sub_test_df))
      all_sex_auc[j,(sex*3+1):(sex*3+3)] <- as.numeric(ci.auc(test_roc))
    }
  }
} else {
  all_sex_auc <- NULL
}

return_sex_split <- list("all_sex_auc" = all_sex_auc)

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

return_ethnic_groups <- list("names_to_keep" = colnames(ethnic_scores), "brit_stats" = brit_stats, "euro_stats" = euro_stats, "african_stats" = african_stats, "asian_stats" = asian_stats)


##############################################
########### SIBLING GROUPS ###################
##############################################

#have to look at sibling pairs compared to non-sibling pairs

sibs <- read.table("../get_covars/covar_data/sibs", stringsAsFactors=F, header=T)

use_pheno <- unname(apply(brit_pheno[,1:4], 1, function(x) 1 %in% x)*1)
df <- data.frame(pheno = use_pheno, brit_covars)

conc_ci <- function(conc, ss, cases){
  #see: The meaning and use of the area under a receiver operating characteristic (ROC) curve
  q1 <- conc/(2-conc)
  q2 <- (2*conc^2)/(1+conc)
  var <- ((conc*(1-conc) + (cases-1)*(q1-conc^2))+(ss-cases-1)*(q2-conc^2))/(cases*(ss-cases))
  c_logit <- log(conc/(1-conc))
  c_var_logit <- var/((conc*(1-conc))^2)
  ci_hi <- exp(c_logit+1.96*c_var_logit)/(1+exp(c_logit+1.96*c_var_logit))
  ci_lo <- exp(c_logit-1.96*c_var_logit)/(1+exp(c_logit-1.96*c_var_logit))
  return(matrix(c(ci_hi, ci_lo)))
}

all_sib_conc <- matrix(0, nrow = ncol(brit_scores), ncol = 6)
for(i in 1:ncol(brit_scores)){
  print(i)
  df$score <- brit_scores[,i]

  pheno_sibs <- sibs[sibs[,1] %in% df$eid[df$pheno == 1] | sibs[,2] %in% df$eid[df$pheno == 1],]
  pos_guess <- rep(NA, nrow(pheno_sibs))
  for(j in 1:nrow(pheno_sibs)){
    sub_df <- df[df$eid %in% pheno_sibs[j,],]
    if(nrow(sub_df) == 2){
      if(sum(sub_df$pheno) == 1){
        pos_guess[j] <- (sub_df$score[sub_df$pheno == 1] > sub_df$score[sub_df$pheno == 0])*1
      }
    }
  }

  pos_guess <- pos_guess[!is.na(pos_guess)]
  sib_conc <- sum(pos_guess)/length(pos_guess)
  if(sib_conc == 1){
    sib_conc_ci <- c(1,1)
  } else {
    sib_conc_ci <- as.numeric(conc_ci(sib_conc, length(pos_guess)*2, length(pos_guess)))
  }

  pos_guess <- rep(NA, nrow(pheno_sibs)*10)
  for(j in 1:length(pos_guess)){
    pos_guess[j] <- (df$score[sample(which(df$pheno == 1), 1)] > df$score[sample(which(df$pheno == 0), 1)])*1
  }

  nonsib_conc <- sum(pos_guess)/length(pos_guess)
  nonsib_conc_ci <- as.numeric(conc_ci(nonsib_conc, length(pos_guess)*2, length(pos_guess)))

  all_sib_conc[i,] <- c(sib_conc_ci[1], sib_conc, sib_conc_ci[2], nonsib_conc_ci[1], nonsib_conc, nonsib_conc_ci[2])
}



rownames(all_sib_conc) <- colnames(brit_scores)
colnames(all_sib_conc) <- c("sib_conc_ci_lo", "sib_conc", "sib_conc_ci_hi", "nonsib_conc_ci_lo", "nonsib_conc", "nonsib_conc_ci_hi")

return_sibling_groups <- list("all_sib_conc" = all_sib_conc)







###########################

extra_info <- list(score_names = colnames(brit_scores))
save_obj <- list("score_sizes" = return_score_sizes, "pheno_defs" = return_pheno_defs, "sex_split" = return_sex_split, "ethnic" = return_ethnic_groups, "sibs" = return_sibling_groups, "extra" = extra_info)
saveRDS(save_obj, paste0("per_score_results/", tolower(author), ".res.RDS"))
