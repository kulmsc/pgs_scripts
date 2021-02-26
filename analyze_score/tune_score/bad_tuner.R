library(survival)
library(pROC)
library(epitools)

#author = commandArgs(trailingOnly=TRUE)
author <- "Nikpay"
#author <- "Christophersen"
phen_method <- "icd_selfrep"
subrate_style <- "slow"
train_frac <- 0.6

input_folds <- 3
input_repeats <- 3


#Read in the PGSs and sort down to training
all_scores <-  readRDS(paste0("../../do_score/final_scores/all_score.", tolower(author), ".RDS"))
all_eid <- read.table("~/athena/doc_score/do_score/all_specs/for_eid.fam", stringsAsFactors=F)
train_eid <- read.table(paste0("~/athena/doc_score/qc/cv_files/train_eid.", train_frac, ".txt"), stringsAsFactors=F)
all_scores <- all_scores[all_eid[,1] %in% train_eid[,1],]
eid <- all_eid[all_eid[,1] %in% train_eid[,1],1]

#Normalize the scores
scores <- all_scores[,grepl(tolower(author), colnames(all_scores))]
scores <- scores[,apply(scores, 2, function(x) length(unique(x)) > 3)]
scores <- apply(scores, 2, function(x) (x-mean(x)) / (max(abs((x-mean(x)))) * 0.01) )


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
#exit()

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
dates <- dates[pheno_eids %in% eid]
pheno <- pheno[pheno_eids %in% eid]
pheno_eids <- pheno_eids[pheno_eids %in% eid]
dates <- dates[order(pheno_eids)[rank(eid)]]
pheno <- pheno[order(pheno_eids)[rank(eid)]]
#Note that pheno_eids should not be used after this point, only eid


#Read in the base covars
covars <- readRDS("../get_covars/base_covars.RDS")
covars <- covars[covars[,1] %in% eid,]
covars <- covars[order(covars[,1])[rank(eid)],]

#Set up survival analysis data frame
#Artifically decide start date is 1999, that way all are even, if date is prior then remove it
#The current maximum date possible is 31 May 2020
death <- read.table("~/athena/ukbiobank/hesin/death.txt", stringsAsFactors=F, header = T)
death[,5] <- unlist(lapply(death[,5], function(x) paste0(strsplit(x, "/")[[1]][3], "-", strsplit(x, "/")[[1]][2], "-", strsplit(x, "/")[[1]][1])))
death <- death[!duplicated(death[,1]),]
death <- death[death[,1] %in% eid,]
add_on <- death[1,]
add_on[5] <- ""
add_eid <- eid[!(eid %in% death[,1])]
add_on <- add_on[rep(1, length(add_eid)),]
add_on$eid <- add_eid
death <- rbind(death, add_on)
death <- death[order(death[,1])[rank(eid)],]

censor <- read.table("../get_covars/covar_data/censor_covars", stringsAsFactors=F, header = T, sep = ",")
censor <- censor[censor[,1] %in% eid,]
censor <- censor[order(censor[,1])[rank(eid)],]


start_date <- rep("1999-01-01", nrow(scores))
end_date <- rep("2020-05-31", nrow(scores))
end_date[censor[,2] != ""] <- censor[censor[,2] != "", 2]
end_date[death[,5] != ""] <- death[death[,5] != "", 5]
end_date[!is.na(dates)] <- dates[!is.na(dates)]

#assigned it to a variable so we can more easily turn off this option
remove_eids <- eid[as.Date(dates) < as.Date("1999-01-01") & !is.na(dates)]

is_death_date <- rep(0, nrow(scores))
is_death_date[death[,5] != ""] <- 1
surv_df <- data.frame(time = as.numeric(as.Date(end_date) - as.Date(start_date)), pheno, is_death_date, covars)
df <- data.frame(pheno, covars[,-1])

#need to remove people where the diagnosis occured before the percieved start date

#need to remove people that had diagnosis before accepted start of the study
scores <- scores[surv_df$time > 0,]
df <- df[surv_df$time > 0,]
covars <- covars[surv_df$time > 0,]
pheno <- pheno[surv_df$time > 0]
surv_df <- surv_df[surv_df$time > 0,]

#author_defs
author_defs <- read.table("../descript_defs/author_defs", stringsAsFactors=F, header=T)
subset_sex <- author_defs$sex[author_defs$author == author]
if(subset_sex == "F"){
  pheno <- pheno[df$sex == 0]
  covars <- covars[df$sex == 0,]
  scores <- scores[df$sex == 0,]
  eid <- eid[df$sex == 0]
  df <- df[df$sex == 0,]
  surv_df <- surv_df[surv_df$sex == 0,]
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
}else if(subset_sex == "M"){
  pheno <- pheno[df$sex == 1]
  covars <- covars[df$sex == 1,]
  scores <- scores[df$sex == 1,]
  eid <- eid[df$sex == 1]
  df <- df[df$sex == 1,]
  surv_df <- surv_df[surv_df$sex == 1,]
  base_covars <- c("age + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
} else {
  base_covars <- c("age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
}

#Set up Fine and Gray
print("finegray")
event_type <- rep("censor", nrow(surv_df))
event_type[surv_df$pheno == 1] <- "diagnosis"
event_type[surv_df$is_death_date == 1] <- "death"
surv_df$event_type <- as.factor(event_type)
fg_diag <- finegray(Surv(time, event_type) ~ ., data = cbind(surv_df[,-which(colnames(surv_df) %in% c("is_death_date", "pheno"))], scores), etype="diagnosis")
fg_death <- finegray(Surv(time, event_type) ~ ., data = cbind(surv_df[,-which(colnames(surv_df) %in% c("is_death_date", "pheno"))], scores), etype="death")
fg_scores_diag <- fg_diag[,grepl(tolower(author), colnames(fg_diag))]
fg_scores_death <- fg_death[,grepl(tolower(author), colnames(fg_death))]
fg_diag <- fg_diag[,!grepl("ss", colnames(fg_diag))]
fg_death <- fg_death[,!grepl("ss", colnames(fg_death))]


#Need to set up repeated cross validation
size_group <- floor(nrow(covars)/input_folds)
start_spots <- floor(seq(1, size_group, length.out=input_repeats+1))
repeat_list <- list()
for(i in 1:input_repeats){
  folds_list <- list()
  for(j in 1:input_folds){
    test_group <- eid[(start_spots[i]+((j-1)*size_group)):(start_spots[i]+(j*size_group))]
    train_group <- eid[!(eid %in% test_group)]
    train_group <- train_group[!(train_group %in% remove_eids)]
    test_group <- test_group[!(test_group %in% remove_eids)]
    folds_list[[j]] <- list("train" = train_group, "test" = test_group)
  }
  repeat_list[[i]] <- folds_list
}


overall_counter <- 1
all_conc_holder <- list()
all_survfit_holder <- list()
all_auc_holder <- list()
all_or_holder <- list()
all_base_holder <- replicate(4, list())

for(nrepeat in 1:input_repeats){
  for(nfold in 1:input_folds){
    print(paste("nfold", nfold))
    print(paste("nrepeat", nrepeat))

    print("survival")
    # Set up the survival analysis data
    fg_diag_train <- fg_diag[fg_diag$eid %in% repeat_list[[nrepeat]][[nfold]][["train"]],]
    fg_diag_test <- fg_diag[fg_diag$eid %in% repeat_list[[nrepeat]][[nfold]][["test"]],]

    surv_df_train <- surv_df[surv_df$eid %in% repeat_list[[nrepeat]][[nfold]][["train"]],]
    surv_df_test <- surv_df[surv_df$eid %in% repeat_list[[nrepeat]][[nfold]][["test"]],]

    
    ###########################################################
    #                    SURVIVAL ANALYSES                    #
    ###########################################################
    
    base_mod <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~", base_covars)), data = fg_diag_train, weight = fgwt)
    
    base_conc <- survConcordance(Surv(fgstart, fgstop, fgstatus) ~ predict(base_mod, fg_diag_test), data = fg_diag_test, weight = fgwt)

    base_survfit <- survfit(base_mod, fg_diag_test, se.fit = F)
    base_full_cumhaz <- base_survfit$cumhaz[,!duplicated(fg_diag_test$eid)]
    final_cumhaz <- base_full_cumhaz[nrow(base_survfit$cumhaz),]
    
    group_factor <- rep(1, length(final_cumhaz))
    group_factor[final_cumhaz < quantile(final_cumhaz, 0.2)] <- 0
    group_factor[final_cumhaz > quantile(final_cumhaz, 0.8)] <- 2
    
    base_cumhaz <- data.frame(time = base_survfit$time,
                                    mean_lo = apply(base_full_cumhaz[,group_factor==0], 1, mean),
                                    mean_mid = apply(base_full_cumhaz[,group_factor==1], 1, mean),
                                    mean_hi = apply(base_full_cumhaz[,group_factor==2], 1, mean),
                                    sd_lo = apply(base_full_cumhaz[,group_factor==0], 1, sd),
                                    sd_mid = apply(base_full_cumhaz[,group_factor==1], 1, sd),
                                    sd_hi = apply(base_full_cumhaz[,group_factor==2], 1, sd))
    base_cumhaz <- tail(base_cumhaz[!duplicated(base_cumhaz$mean_mid),],1)
    rm(base_survfit)
    rm(base_full_cumhaz)
    
    all_conc <- matrix(0, nrow = ncol(scores), ncol = 2)
    all_subrates <- list()
    
    for(i in 1:nrow(all_conc)){
      print(i)
      fg_diag_train$score <- fg_scores_diag[fg_diag$eid %in% repeat_list[[nrepeat]][[nfold]][["train"]],i]
      fg_diag_test$score <- fg_scores_diag[fg_diag$eid %in% repeat_list[[nrepeat]][[nfold]][["test"]],i]
      surv_df_train$score <- scores[surv_df$eid %in% repeat_list[[nrepeat]][[nfold]][["train"]],i]
      surv_df_test$score <- scores[surv_df$eid %in% repeat_list[[nrepeat]][[nfold]][["test"]],i]

      #Make the model
      score_mod <- coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~", base_covars, "+ score")), data = fg_diag_train, weight = fgwt)
      
      #CONCORDANCE ###################################
      score_conc_obj <- survConcordance(Surv(fgstart, fgstop, fgstatus) ~ predict(score_mod, fg_diag_test), data = fg_diag_test, weight = fgwt)

      all_conc[i,1] <- score_conc_obj$conc
      all_conc[i,2] <- score_conc_obj$std.err
      
     
      exit()
      #SURVFIT ##########################################
      #SLOW
      if(subrate_style == "slow"){

        H0 <- basehaz(score_mod, centered = FALSE)
        coef <- score_mod$coefficients
        good_order <- colnames(surv_df_test)[colnames(surv_df_test) %in% names(coef)]
        coef <- coef[order(names(coef))[rank(good_order)]]
        prop_haz <- rowSums(t(t(surv_df_test[,colnames(surv_df_test) %in% names(coef)]) * coef))
        H <- H0$hazard[nrow(H0)] * exp(prop_haz)








        split_group <- split(unique(fg_diag_test$eid), rep(1:ceiling(length(unique(fg_diag_test$eid))/100), each = 100))

        get_smaller <- function(smaller_test){
          score_full_cumhaz <- survfit(score_mod, smaller_test, se.fit = F)
          later_time <- score_full_cumhaz$time
          score_full_cumhaz <- score_full_cumhaz$cumhaz
          score_full_cumhaz <- score_full_cumhaz[,!duplicated(smaller_test$eid)]
          final_cumhaz <- score_full_cumhaz[nrow(score_full_cumhaz),]
        
          return(final_cumhaz)
        }

        get_end <- function(score_full_cumhaz){
          group_factor <- rep(1, length(score_full_cumhaz))
          group_factor[score_full_cumhaz < quantile(score_full_cumhaz, 0.2)] <- 0
          group_factor[score_full_cumhaz > quantile(score_full_cumhaz, 0.8)] <- 2
        
          pred_cumhaz_score <- data.frame(mean_lo = mean(score_full_cumhaz[group_factor==0]),
                                        mean_mid = mean(score_full_cumhaz[group_factor==1]),
                                        mean_hi = mean(score_full_cumhaz[group_factor==2]),
                                        sd_lo = sd(score_full_cumhaz[group_factor==0]),
                                        sd_mid = sd(score_full_cumhaz[group_factor==1]),
                                        sd_hi = sd(score_full_cumhaz[group_factor==2]))
          return(pred_cumhaz_score)
        }

        temp_list <- list()
        for(ii in 1:length(split_group)){
          temp_list[[ii]] <- get_smaller(fg_diag_test[fg_diag_test$eid %in% split_group[[ii]],])
        }
        pred_cumhaz_score <- get_end(unlist(temp_list))

      } else if(subrate_style == "fast"){ 
        #FAST
        group_list <- list()
        for(k in 1:length(names(score_mod$coef))){
         group_list[[k]] <-  as.numeric(quantile(surv_df_train[[names(score_mod$coef)[k]]], c(0.1, 0.5, 0.9)))
         if(sign(score_mod$coef[k] == -1)){
          group_list[[k]] <- rev(group_list[[k]])
         }
        }
        mean_df <- data.frame(do.call("cbind", group_list))
        colnames(mean_df) <- names(score_mod$coef)
        if(sign(score_mod$coef[2]) == -1){
         mean_df$sex <- c(0.9, 0.5, 0.1)
        } else {
         mean_df$sex <- c(0.1, 0.5, 0.9)
        }
        
        score_pred <- survfit(score_mod, newdata = mean_df)
        
        pred_cumhaz_score <- data.frame(score_pred$time, score_pred$cumhaz, score_pred$std.err)
        colnames(pred_cumhaz_score) <- c("time", "mean_lo", "mean_mid", "mean_hi", "sd_lo", "sd_mid", "sd_hi")
        pred_cumhaz_score <- pred_cumhaz_score[!duplicated(pred_cumhaz_score$mean_mid),]
      }
      
      
      all_subrates[[i]] <- pred_cumhaz_score
      
    }





    # Set up the normal model data
    print("normal")

    df_train <- df[covars[,1] %in% repeat_list[[nrepeat]][[nfold]][["train"]],]
    df_test <- df[covars[,1] %in% repeat_list[[nrepeat]][[nfold]][["test"]],]
    scores_train <- scores[covars[,1] %in% repeat_list[[nrepeat]][[nfold]][["train"]],]
    scores_test <- scores[covars[,1] %in% repeat_list[[nrepeat]][[nfold]][["test"]],]
    pheno_test <- pheno[covars[,1] %in% repeat_list[[nrepeat]][[nfold]][["test"]]]
  
    ###########################################################
    #                    NORMAL MODELS                        #
    ###########################################################

    base_mod <- glm(pheno ~ ., data = df_train, family = "binomial")
    base_pred <- predict(base_mod, df_test)
    base_roc <- roc(pheno_test ~ base_pred, quiet=T)
    base_auc <- as.numeric(ci.auc(base_roc))

    base_group <- rep(1, length(base_pred))
    base_group[base_pred < quantile(base_pred, 0.2)] <- 0
    base_group[base_pred > quantile(base_pred, 0.8)] <- 2
    base_odds_table <- matrix(c(sum(df_test$pheno == 1 & base_group == 2), sum(df_test$pheno == 0 & base_group == 2),
                         sum(df_test$pheno == 1 & base_group == 0), sum(df_test$pheno == 0 & base_group == 0)), nrow = 2)
    base_odds_ratio <- oddsratio.wald(base_odds_table)
    base_odds_ratio <- base_odds_ratio$measure[2,c(2,1,3)]
  
    all_auc <- matrix(0, nrow = ncol(scores), ncol = 3)
    all_odds_ratio <- matrix(0, nrow = ncol(scores), ncol = 3)
    for(i in 1:ncol(scores)){
      df_train$score <- scores_train[,i]
      df_test$score <- scores_test[,i]
      #Make the model
      score_mod <- glm(pheno ~ ., data = df_train, family = "binomial")
      score_pred <- predict(score_mod, df_test)
      
      # AUC #####################################
      score_roc <- roc(pheno_test ~ score_pred, quiet=T)
      all_auc[i,] <- as.numeric(ci.auc(score_roc))
  
      # ODDS RATIO #############################
      score_group <- rep(1, length(score_pred))
      score_group[score_pred < quantile(score_pred, 0.2)] <- 0
      score_group[score_pred > quantile(score_pred, 0.8)] <- 2
  
      score_odds_table <- matrix(c(sum(df_test$pheno == 1 & score_group == 2), sum(df_test$pheno == 0 & score_group == 2),
                         sum(df_test$pheno == 1 & score_group == 0), sum(df_test$pheno == 0 & score_group == 0)), nrow = 2)
      score_odds_ratio <- oddsratio.wald(score_odds_table)
      all_odds_ratio[i,] <- score_odds_ratio$measure[2,c(2,1,3)]
    }
    
    all_conc_holder[[overall_counter]] <- all_conc
    all_survfit_holder[[overall_counter]] <- all_subrates
    all_auc_holder[[overall_counter]] <- all_auc
    all_or_holder[[overall_counter]] <- all_odds_ratio

    all_base_holder[[1]][[overall_counter]] <- base_conc
    all_base_holder[[2]][[overall_counter]] <- base_cumhaz
    all_base_holder[[3]][[overall_counter]] <- base_auc
    all_base_holder[[4]][[overall_counter]] <- base_odds_ratio
    overall_counter <- overall_counter + 1
  }
}

misc_info <- list("phen_method" = phen_method, "subrate_style" = subrate_style, "train_frac" = train_frac)
final_obj <- list("conc" = all_conc_holder, "survfit" = all_survfit_holder, "auc" = all_auc_holder, "or" = all_or_holder, "base" = all_base_holder, "score_names" = colnames(scores), "misc" = misc_info, "overall_df" = df, "surv_df" = surv_df)
saveRDS(surv_df, paste0("tune_results/", author, "_survdf.RDS"))
saveRDS(final_obj, paste0("tune_results/", author, "_res.RDS"))


