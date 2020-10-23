library(survival)
library(pROC)
library(epitools)

author <- "Bentham"
phen_method <- "all"

#Read in the PGSs and sort down to training
all_scores <- list()
all_files <- list.files("~/athena/doc_score/do_score/final_scores", "RDS")
for(i in 1:length(all_files)){
  all_scores[[i]] <- readRDS(paste0("~/athena/doc_score/do_score/final_scores/", all_files[i]))
}
all_scores <- do.call("cbind", all_scores)

pheno_eids <- read.table("../construct_defs/eid.csv", header = T)
all_eid <- read.table("~/athena/doc_score/do_score/all_specs/for_eid.fam", stringsAsFactors=F)
train_eid <- read.table("~/athena/doc_score/qc/cv_files/train_eid.0.2.txt", stringsAsFactors=F)
test_eid <- read.table("~/athena/doc_score/qc/cv_files/test_eid.0.8.txt", stringsAsFactors=F)
pheno_eids <- read.table("../construct_defs/eid.csv", header = T)

#Normalize the scores
best_name <- read.table(paste0("../tune_score/tune_results/", tolower(author), ".best.ss"), stringsAsFactors=F)
best_score <- all_scores[,colnames(all_scores) == best_name[1,1]]
best_score <- (best_score-mean(best_score)) / (max(abs((best_score-mean(best_score)))) * 0.01)
train_pgs <- best_score[all_eid[,1] %in% train_eid[,1] & all_eid[,1] %in% pheno_eids[,1]]
test_pgs <- best_score[all_eid[,1] %in% test_eid[,1] & all_eid[,1] %in% pheno_eids[,1]]
rm(all_scores)


#Read in the phenotypes, order is: cancer sr, noncancer sr, icd9, icd10, oper, meds
#selfasses date: 2018-11-22; hesin date: 21/01/2000
pheno <- read.table(paste0("../construct_defs/pheno_defs/diag.", tolower(author), ".txt.gz"), stringsAsFactors=F)
dates <- read.table(paste0("../construct_defs/pheno_defs/time.", tolower(author), ".txt.gz"), stringsAsFactors=F)
for(i in 1:ncol(dates)){
  if(i %in% c(1,2,6)){
    dates[dates[,i] == "__________", i] <- "2020-12-31"
    dates[,i] <- as.Date(dates[,i], "%Y-%m-%d")
  } else {
    dates[dates[,i] == "__________", i] <- "31/12/2020"
    dates[,i] <- as.Date(dates[,i], "%d/%m/%Y")
  }
}

if(phen_method == "icd"){
  pheno <- pheno[,3:4]
  dates <- dates[,3:4]
} else if(phen_method == "selfrep"){
  pheno <- pheno[,1:2]
  dates <- dates[,1:2]
} else if(phen_method == "icd_selfrep"){
  pheno <- pheno[,1:4]
  dates <- dates[,1:4]
} else if(phen_method == "all" | phen_method == "double"){
  print("doing nothing")
}

dates <- apply(dates, 1, min)
dates[dates == as.Date("2020-12-31")] <- NA
if(phen_method == "double"){
  pheno <- rowSums(pheno)
  pheno[pheno == 1] <- 0
  pheno[pheno > 1] <- 1

  dates[pheno == 0] <- NA
} else {
  pheno <- rowSums(pheno)
  pheno[pheno > 1] <- 1
}


#Read in the eids used that are the same order as the pheno and dates, then subset the pheno and dates accordingly
pheno_eids <- read.table("../construct_defs/eid.csv", header = T)
pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
pheno_eids <- pheno_eids[-length(pheno_eids)]
dates <- dates[pheno_eids %in% all_eid[,1]]
pheno <- pheno[pheno_eids %in% all_eid[,1]]
pheno_eids <- pheno_eids[pheno_eids %in% all_eid[,1]]
dates <- dates[order(pheno_eids)[rank(all_eid[,1])]]
pheno <- pheno[order(pheno_eids)[rank(all_eid[,1])]]

#Read in the base covars
covars <- readRDS("../get_covars/base_covars.RDS")
covars <- covars[covars[,1] %in% all_eid[,1],]
covars <- covars[order(covars[,1])[rank(all_eid[,1])],]

#Set up survival analysis data frame
#Artifically decide start date is 1999, that way all are even, if date is prior then remove it
#The current maximum date possible is 31 May 2020
death <- read.table("~/athena/ukbiobank/hesin/death.txt", stringsAsFactors=F, header = T)
death[,5] <- unlist(lapply(death[,5], function(x) paste0(strsplit(x, "/")[[1]][3], "-", strsplit(x, "/")[[1]][2], "-", strsplit(x, "/")[[1]][1])))
death <- death[!duplicated(death[,1]),]
death <- death[death[,1] %in% all_eid[,1],]
add_on <- death[1,]
add_on[5] <- ""
add_eid <- all_eid[!(all_eid[,1] %in% death[,1]),1]
add_on <- add_on[rep(1, length(add_eid)),]
add_on$eid <- add_eid
death <- rbind(death, add_on)
death <- death[order(death[,1])[rank(all_eid[,1])],]

censor <- read.table("../get_covars/covar_data/censor_covars", stringsAsFactors=F, header = T, sep = ",")
censor <- censor[censor[,1] %in% all_eid[,1],]
if(nrow(all_eid) > nrow(censor)){
  add_on_eid <- all_eid[!all_eid[,1] %in% censor[,1],1]
  add_on <- censor[rep(1, length(add_on_eid)),]
  add_on[,1] <- add_on_eid
  censor <- rbind(censor, add_on)
}
censor <- censor[order(censor[,1])[rank(all_eid[,1])],]

start_date <- rep("1999-01-01", length(best_score))
end_date <- rep("2020-05-31", length(best_score))
end_date[censor[,2] != ""] <- censor[censor[,2] != "", 2]
end_date[death[,5] != ""] <- death[death[,5] != "", 5]

is_death_date <- rep(0, length(best_score))
is_death_date[death[,5] != ""] <- 1
surv_df <- data.frame(time = as.numeric(as.Date(end_date) - as.Date(start_date)), pheno, is_death_date, covars, best_score)
surv_df_train <- surv_df[surv_df$eid %in% train_eid[,1],]
surv_df_test <- surv_df[surv_df$eid %in% test_eid[,1],]


#Set up Fine and Gray
print("finegray")

event_type <- rep("censor", nrow(surv_df_train))
event_type[surv_df_train$pheno == 1] <- "diagnosis"
event_type[surv_df_train$is_death_date == 1] <- "death"
surv_df_train$event_type <- as.factor(event_type)
fg_diag_train <- finegray(Surv(time, event_type) ~ ., data = surv_df_train[,-which(colnames(surv_df_train) %in% c("is_death_date", "eid", "pheno"))], etype="diagnosis")
fg_death_train <- finegray(Surv(time, event_type) ~ ., data = surv_df_train[,-which(colnames(surv_df_train) %in% c("is_death_date", "eid", "pheno"))], etype="death")

event_type <- rep("censor", nrow(surv_df_test))
event_type[surv_df_test$pheno == 1] <- "diagnosis"
event_type[surv_df_test$is_death_date == 1] <- "death"
surv_df_test$event_type <- as.factor(event_type)
fg_diag_test <- finegray(Surv(time, event_type) ~ ., data = surv_df_test[,-which(colnames(surv_df_test) %in% c("is_death_date", "eid", "pheno"))], etype="diagnosis")
fg_death_test <- finegray(Surv(time, event_type) ~ ., data = surv_df_test[,-which(colnames(surv_df_test) %in% c("is_death_date", "eid", "pheno"))], etype="death")


#Get Conc
print("CONC")
base_mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = fg_diag_train, weight = fgwt)
score_mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + best_score, data = fg_diag_train, weight = fgwt)

base_conc_obj <- survConcordance(Surv(fgstart, fgstop, fgstatus) ~ predict(base_mod, fg_diag_test), data = fg_diag_test, weight = fgwt)
score_conc_obj <- survConcordance(Surv(fgstart, fgstop, fgstatus) ~ predict(base_mod, fg_diag_test), data = fg_diag_test, weight = fgwt)


#Get the adjusted hazard rates
print("SUBRATES")
#Set up the basic risk groups based on bottom two and top two deciles
train_score_factor <- rep(0, length(train_pgs))
train_score_factor[train_pgs > quantile(train_pgs, 0.2) & train_pgs <= quantile(train_pgs, 0.8)] <- 1
train_score_factor[train_pgs > quantile(train_pgs, 0.8)] <- 2

test_score_factor <- rep(0, length(test_pgs))
test_score_factor[test_pgs > quantile(test_pgs, 0.2) & test_pgs <= quantile(test_pgs,0.8)] <- 1
test_score_factor[test_pgs > quantile(test_pgs, 0.8)] <- 2

fg_pgs <- fg_diag_train$best_score
fg_score_factor <- rep(0, length(fg_pgs))
fg_score_factor <- rep(0, length(fg_pgs))
fg_score_factor[fg_pgs > quantile(fg_pgs, 0.2) & fg_pgs <= quantile(fg_pgs,0.8)] <- 1
fg_score_factor[fg_pgs > quantile(fg_pgs, 0.8)] <- 2

##
#cannote use predict, as it outputs hazard ratios with no time correspondence
#could try to transfer to using accelerated rate time models, but also very complex likely
#can try to do matched case/control to implicitly adjust
#keep reading the vignette

exit()

#does not work
base_mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = fg_diag_train, weight = fgwt)
#fg_survfit <- survfit(base_mod, newdata = fg_diag_test, data = fg_diag_train, weight = fgwt)
fg_survfit <- survfit(base_mod, fg_diag_test[1:10,])
fg_survfit2 <- survfit(base_mod, fg_diag_test[1:10,], weight=fgwt)
fg_survfit3 <- survfit(base_mod, data=fg_diag_train, new_data=fg_diag_test[1:10,], weight=fgwt)
fg_survfit4 <- survfit(base_mod, surv_df_test[1:10,])
#should try manually changing the coefs of the simple_base_mod to those of base_mod

#setting up simple models
simple_base_mod <- coxph(Surv(time, pheno) ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = surv_df_train)
#base_survfit <- survfit(simple_base_mod, newdata = surv_df_test, data = surv_df_train)
base_survfit <- survfit(simple_base_mod, surv_df_test[1:10,])
base_survfit2 <- survfit(simple_base_mod, fg_diag_test[1:10,])
base_survfit3 <- survfit(simple_base_mod, fg_diag_test[1:10,], weight = fgwt)
base_survfit4 <- survfit(simple_base_mod, data=surv_df_train, newdata=fg_diag_test[1:10,], weight = fgwt)
simple_base_mod$coef <- base_mod$coef
base_survfit5 <- survfit(simple_base_mod, surv_df_test[1:10,])

#setting up multievent
surv_df_train$multievent <- 0
surv_df_train$multievent[surv_df_train$pheno == 1] <- 1
surv_df_train$multievent[surv_df_train$is_death_date == 1] <- 2
surv_df_test$multievent <- 0
surv_df_test$multievent[surv_df_test$pheno == 1] <- 1
surv_df_test$multievent[surv_df_test$is_death_date == 1] <- 2
surv_df_train$multievent <- as.factor(surv_df_train$multievent)
surv_df_test$multievent <- as.factor(surv_df_test$multievent)

#takes way too long
multi_mod <- coxph(Surv(time, multievent) ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = surv_df_train, id = eid)
multi_survfit <- survfit(multi_mod, newdata = surv_df_test[1:1000,], data = surv_df_train)

#set up stratification
surv_df_train$score_group <- as.factor(train_score_factor+1)
surv_df_test$score_group <- as.factor(test_score_factor+1)

#do the simple stratified model
strata_base_mod <- coxph(Surv(time, pheno) ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + strata(score_group), data = surv_df_train)
easy_fit <- survfit(simple_base_mod, newdata = surv_df_test)

#simple stratified but with mutimodel - does not work
strata_multi_mod <- coxph(Surv(time, multievent) ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + strata(score_group), data = surv_df_train, id = eid)
easy_fit <- survfit(strata_multi_mod)

#simple stratified but with fg
fg_diag_train$score_group <- as.factor(fg_score_factor)
strata_fg_mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~  age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + strata(score_group), data = fg_diag_train, weight = fgwt)
easy_fit <- survfit(strata_fg_mod)

#carry out vignette approach, does not work
tdata <- data.frame(surv_df_train[1:10000,5:16])
tdata3 <- tdata[rep(1:nrow(tdata), 3),] #three copies
tdata3$score_group <- factor(rep(1:3, each=nrow(tdata)), levels = levels(surv_df_train$score_group))
sfit4a <- survexp(~score_group, data=tdata3, ratetable=simple_base_mod)

#new plan, do multistate then fill in with dummy as done in the vignette, if possible expand dummy to include more of the people in the dataset and then average
multi_mod <- coxph(Surv(time, multievent) ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + best_score, data = surv_df_train, id = eid)
avg_row <- apply(surv_df_test[,5:16], 2, mean)
dummy <- as.data.frame(t(data.frame(avg_row))[rep(1,3),])
dummy$best_score <- as.numeric(quantile(surv_df_test$best_score, c(0.1, 0.5, 0.9)))
ans=survfit(multi_mod, dummy)





###

  adjrates <- data.frame(time = rep(score_fit$time, 3), n = rep(score_fit$n.risk, 3), cumhaz = as.numeric(score_fit$cumhaz),
     stderr = as.numeric(score_fit$std.err), group = rep(1:3, each = length(score_fit$time)))
  subrates <- adjrates[adjrates$time == max(adjrates$time),]
  all_subrates[[i]] <- subrates


#AUC
print("AUC")
df <- data.frame(pheno, covars[,-1])
base_mod <- glm(pheno ~ ., data = df, family = "binomial")
base_roc <- roc(pheno ~ predict(base_mod, df))
base_auc <- as.numeric(ci.auc(base_roc))
all_auc <- matrix(0, nrow = ncol(scores), ncol = 3)
for(i in 1:ncol(scores)){
  df$score <- scores[,i]
  score_mod <- glm(pheno ~ ., data = df, family = "binomial")
  score_roc <- roc(pheno ~ predict(score_mod, df))
  all_auc[i,] <- as.numeric(ci.auc(score_roc))
}

#ODDS RATIOS
print("ODDS RATIO")
all_odds_ratio <- matrix(0, nrow = ncol(scores), ncol = 3)
for(i in 1:ncol(scores)){
  df$group <- 1
  df$group[scores[,i] < quantile(scores[,i], 0.2)] <- 0
  df$group[scores[,i] > quantile(scores[,i], 0.8)] <- 2
  odds_table <- matrix(c(sum(df$pheno == 1 & df$group == 2), sum(df$pheno == 0 & df$group == 2),
                       sum(df$pheno == 1 & df$group == 0), sum(df$pheno == 0 & df$group == 0)), nrow = 2)
  odds_ratio <- oddsratio.wald(odds_table)
  all_odds_ratio[i,] <- odds_ratio$measure[2,c(2,1,3)]
}

#PICK THE BEST
if(choose_stat == "conc"){
  best_ss <- colnames(scores)[which.max(all_conc[,1])]
} else if(choose_stat == "auc"){
  best_ss <- colnames(scores)[which.max(all_auc[,2])]
} else if(choose_stat == "or"){
  best_ss <- colnames(scores)[which.max(all_odds_ratio[,2])]
} else if(choose_stat == "cumhaz"){
  best_ss <- colnames(scores)[which.max(unlist(lapply(all_subrates, function(x) x[3,3]/x[1,3])))]
}

save_list <- list(eid = eid, author = author, phen_method = phen_method, pheno = pheno, 
                 dates = dates, all_conc = all_conc, all_subrates = all_subrates, all_auc = all_auc,
                 all_odds_ratio = all_odds_ratio, score_names = colnames(scores),
                 base_conc = base_conc, base_auc = base_auc)
saveRDS(save_list, paste0("tune_results/", tolower(author), ".tune.RDS"))

write.table(best_ss, paste0("tune_results/", tolower(author), ".best.ss"), row.names = F, col.names = F, quote = F)
