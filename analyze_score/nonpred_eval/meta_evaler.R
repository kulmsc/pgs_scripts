
#When all of the tuning is done want to go though and check that 
#If method performance is dependent on any factors (sample size, SNPs, etc.)
  #Do for one AUC for each author (the best)
  #Do this once for each method (the best method for each author)

#Want to see if there is a correlation between method parameters and the meta stats
  #so collect all results (for all authors) then run a model

#Also want to compare training to testing results

#Will focus on AUC

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


options(warn = 2)

list_all_auc <- list()
list_all_best_names <- list()
list_all_best_methods <- list()
list_ukb_cc <- list()
list_ukb_snps <- list()
list_ukb_beta_pdiff <- list()
list_ukb_len_pdiff <- list()
list_pval_6 <- list()
list_pval_8 <- list()
list_beta_len <- list()
list_beta_effect <- list()



i <- 1
all_author <- unlist(lapply(strsplit(list.files("../tune_score/tune_results/", "res"), "_"), function(x) x[1]))
all_author <- all_author[all_author != "Xie"]
for(author in all_author){
  all_res <- readRDS(paste0("../tune_score/tune_results/", author, "_res.RDS"))
  score_names <- all_res[["score_names"]]

  #GWAS Stats
  orig_gwas <- read.table(paste0("~/athena/doc_score/raw_ss/", author, "/clean_", tolower(author), ".txt.gz"), stringsAsFactors=F, header=T)
  four_len <- unlist(lapply(splitter(sort(abs(orig_gwas$BETA[orig_gwas$P < 0.05])), 4), function(x) length(x)))
  four_mean <- unlist(lapply(splitter(sort(abs(orig_gwas$BETA[orig_gwas$P < 0.05])), 4), function(x) mean(x)))
  list_beta_len[[i]] <- (four_len[4] - four_len[1])/four_len[1]
  list_beta_effect[[i]] <- (four_mean[4] - four_mean[1])/four_mean[1]
  list_pval_6[[i]] <- sum(orig_gwas$P < 1e-6)
  list_pval_8[[i]] <- sum(orig_gwas$P < 1e-8)

  #read in
  list_all_best_names[[i]] <- read.table(paste0("../tune_score/tune_results/", tolower(author), ".best.ss"), stringsAsFactors=F)
  list_all_best_methods[[i]] <- read.table(paste0("../tune_score/tune_results/", tolower(author), ".methods.ss"), stringsAsFactors=F)

  #get the mean_auc, and put in other info
  mean_auc <- as.data.frame(Reduce("+", all_res[["auc"]])/length(all_res[["auc"]]))
  mean_auc$name <- score_names
  mean_auc$method <- unlist(lapply(strsplit(score_names, ".", fixed=T), function(x) x[3]))
  mean_auc$author <- author

  #get the case control in UKBB
  ukb_pheno <- read.table(paste0("../construct_defs/pheno_defs/diag.", tolower(author) ,".txt.gz"), stringsAsFactors=F)
  list_ukb_cc[[i]] <- t(matrix(as.numeric(table(rowSums(ukb_pheno[,1:4]) > 0))))[rep(1, nrow(mean_auc)),]

  #get the ukb snp lengths
  ukb_snps <- rep(0, nrow(list_all_best_methods[[i]]))
  ukb_beta_pdiff <- rep(0, nrow(list_all_best_methods[[i]]))
  ukb_len_pdiff <- rep(0, nrow(list_all_best_methods[[i]]))
  for(sname in list_all_best_methods[[i]][,1]){
    new_name <- paste0(strsplit(sname, ".", fixed = T)[[1]][3], ".", strsplit(sname, ".", fixed = T)[[1]][2], ".ss.gz")
    get_files <- list.files(paste0("../../mod_sets/", author, "/"), pattern = new_name)
    prs_sets <- list()
    #temp_length <- rep(0, length(get_files))
    for(j in 1:length(get_files)){
      #temp_length[j] <- as.numeric(system(paste0("cat ../../mod_sets/", author, "/", get_files[j], " | wc -l"), intern = TRUE))
      prs_sets[[j]] <- read.table(paste0("../../mod_sets/", author, "/", get_files[j]), stringsAsFactors=F, header=T)
      colnames(prs_sets[[j]]) <- paste0(rep("X", ncol(prs_sets[[j]])), 1:ncol(prs_sets[[j]]))
    }
    prs_sets <- do.call("rbind", prs_sets)

    ukb_snps[which(list_all_best_methods[[i]][,1] == sname)] <- nrow(prs_sets)
    if(length(unique(prs_sets[,7])) > 10){
      temp <- unlist(lapply(splitter(sort(abs(sort(prs_sets[,7]))), 4), function(x) mean(x)))
      ukb_beta_pdiff[which(list_all_best_methods[[i]][,1] == sname)] <- (temp[4] - temp[1])/temp[1]
      temp <- unlist(lapply(splitter(sort(abs(sort(prs_sets[,7]))), 4), function(x) length(x)))
      ukb_len_pdiff[which(list_all_best_methods[[i]][,1] == sname)] <- (temp[4] - temp[1])/temp[1]
    } else {
      ukb_beta_pdiff[which(list_all_best_methods[[i]][,1] == sname)] <- NA
      ukb_len_pdiff[which(list_all_best_methods[[i]][,1] == sname)] <- NA
    }
  }

  names(ukb_snps) <- list_all_best_methods[[i]][,1]
  list_ukb_snps[[i]] <- ukb_snps
  list_ukb_beta_pdiff[[i]] <- ukb_beta_pdiff
  list_ukb_len_pdiff[[i]] <- ukb_len_pdiff


  list_all_auc[[i]] <- mean_auc
  i <- i + 1
}

for(i in 1:length(list_pval_6)){
  list_pval_6[[i]] <- rep(list_pval_6[[i]], length(list_ukb_snps[[i]]))
  list_pval_8[[i]] <- rep(list_pval_8[[i]], length(list_ukb_snps[[i]]))
  list_beta_len[[i]] <- rep(list_beta_len[[i]], length(list_ukb_snps[[i]]))
  list_beta_effect[[i]] <- rep(list_beta_effect[[i]], length(list_ukb_snps[[i]]))
}


#add meta info to the all_auc
all_auc <- do.call("rbind", list_all_auc)
all_auc <- cbind(all_auc, do.call("rbind", list_ukb_cc))

method_best_names <- unlist(list_all_best_methods)
all_auc <- all_auc[all_auc$name %in% method_best_names,]
all_auc <- cbind(all_auc, unlist(list_ukb_snps), unlist(list_ukb_beta_pdiff), unlist(list_ukb_len_pdiff), unlist(list_pval_6), unlist(list_pval_8), unlist(list_beta_len), unlist(list_beta_effect))

colnames(all_auc) <- c("auc_lo", "auc", "auc_hi", "name", "method", "author", "uk_control", "uk_case", "uk_snps", "uk_beta_pdiff", "uk_len_pdiff", "sig6_gwas_snps", "sig8_gwas_snps", "len_pdiff", "effect_pdiff")

meta_info <- read.table("../../raw_ss/meta_stats", stringsAsFactors=F, sep=",", header=T)
add_on <- data.frame(matrix(0, nrow = nrow(all_auc), ncol = 10))
colnames(add_on) <- colnames(meta_info)[-c(1,2)]
for(author in unique(all_auc$author)){
  add_on[all_auc$author == author,] <- meta_info[meta_info$author == author,-c(1,2)]
}
all_auc <- cbind(all_auc, add_on)

all_auc$uk_cc_ratio <- all_auc$uk_case/all_auc$uk_control
all_auc$gwas_cc_ratio <- all_auc$cases/all_auc$controls
all_auc$uk_sampe_size <- all_auc$uk_case + all_auc$uk_control

#get the best names (absolute best, and best name for each method for each author)
abs_best_names <- unlist(lapply(list_all_best_names, function(x) x[3,2]))

for(reg_name in c("uk_snps", "uk_beta_pdiff", "uk_len_pdiff", "sig6_gwas_snps", "sig8_gwas_snps", "len_pdiff", "effect_pdiff", "sampe_size", "snps", "h2", "uk_cc_ratio", "gwas_cc_ratio")){
  all_auc[reg_name] <- (all_auc[reg_name] - min(all_auc[reg_name], na.rm=T))/(max(all_auc[reg_name],na.rm=T) - min(all_auc[reg_name],na.rm=T))
}

#subset the all_auc
all_auc$ancestry <- as.numeric(as.factor(all_auc$ancestry))
all_auc <- all_auc[,-which(colnames(all_auc) %in% c("ldsc_h2", "ldsc_h2_se", "hdl_h2", "hdl_h2_se", "hdl_h2se"))]
abs_best_auc <- all_auc[all_auc$name %in% abs_best_names,]
method_best_auc <- all_auc


exit()
############################################################
################ META INFO REGRESSIONS #####################
############################################################

meta_to_regress <- c("uk_snps", "uk_beta_pdiff", "uk_len_pdiff", "sig6_gwas_snps", "sig8_gwas_snps", "len_pdiff", "effect_pdiff", "sampe_size", "ancestry", "snps", "h2", "uk_cc_ratio", "gwas_cc_ratio")
abs_best_coefs <- matrix(0, nrow = length(meta_to_regress), ncol = 4)
all_poss_coefs <- matrix(0, nrow = length(meta_to_regress), ncol = 4)
method_best_coefs <- rep(list(matrix(0, nrow = length(meta_to_regress), ncol = 4)), length(unique(method_best_auc$method)))

for(i in 1:length(meta_to_regress)){
  best_mod <- lm(as.formula(paste0("auc ~ ", meta_to_regress[i])), data = abs_best_auc)
  abs_best_coefs[i,] <- summary(best_mod)$coefficients[2,]

  all_mod <- lm(as.formula(paste0("auc ~ ", meta_to_regress[i])), data = method_best_auc)
  all_poss_coefs[i,] <- summary(all_mod)$coefficients[2,]

  for(m in unique(method_best_auc$method)){
     method_mod <- lm(as.formula(paste0("auc ~ ", meta_to_regress[i])), data = method_best_auc[method_best_auc$method == m,])
     method_best_coefs[[which(unique(method_best_auc$method) == m)]][i,] <- summary(method_mod)$coefficients[2,]
  }
}


#should add regressions across all for best method, not just abs best auc
meta_to_regress <- meta_to_regress[-which(meta_to_regress == "ancestry")]
umethod <- unique(method_best_auc$method)
method_holder <- data.frame(matrix(0, nrow = length(meta_to_regress)*2, ncol = length(umethod)))
norm_auc <- (method_best_auc$auc - min(method_best_auc$auc))/(max(method_best_auc$auc) - min(method_best_auc$auc))
k <- 1
for(i in 1:length(meta_to_regress)){
  j <- which(colnames(method_best_auc) == meta_to_regress[i])
  stat_val <- (method_best_auc[,j] - min(method_best_auc[,j]))/(max(method_best_auc[,j]) - min(method_best_auc[,j]))
  temp <- method_best_auc[order(norm_auc*stat_val,decreasing=T),]
  top_method <- temp$method[!duplicated(temp$author)][1:10]

  for(m in unique(top_method)){
    method_holder[k,which(umethod == m)] <- sum(top_method == m)
  }

  stat_val <- abs(stat_val - 1)
  temp <- method_best_auc[order(norm_auc*stat_val,decreasing=T),]
  top_method <- temp$method[!duplicated(temp$author)][1:10]

  for(m in unique(top_method)){
    method_holder[k+1,which(umethod == m)] <- sum(top_method == m)
  }

  k <- k + 2
}

colnames(method_holder) <- umethod
method_holder$stat <- paste0(c("", "rev_"), rep(meta_to_regress, each=2))

#want to now determine best method for each type of statistic
#normalize so that each stat is smooshed to 0 to 1, with 1 being best possible
#then multiply auc by the normalized vals and rank


############################################################
################ TRAIN AND TEST COMP #####################
############################################################

ttcomp <- matrix(0, nrow = length(all_author), ncol = 4)

k <- 1
for(author in all_author){
  tune_res <- readRDS(paste0("../tune_score/tune_results/", author, "_res.RDS"))
  test_res <- readRDS(paste0("../test_score/final_stats/", author, "_res.RDS"))
  use_name <- test_res[["score_names"]]
  use_index <- which(tune_res[["score_names"]] == use_name)

  mean_conc <- Reduce("+", tune_res[["conc"]])/length(tune_res[["conc"]])
  mean_auc <- Reduce("+", tune_res[["auc"]])/length(tune_res[["auc"]])
  mean_or <- Reduce("+", tune_res[["or"]])/length(tune_res[["or"]])
  mean_survfit <- matrix(0, nrow = length(tune_res[["survfit"]][[1]]), ncol = 6)
  for(i in 1:length(tune_res[["survfit"]])){
    for(j in 1:length(tune_res[["score_names"]])){
      mean_survfit[j,] <- mean_survfit[j,] + as.numeric(tune_res[["survfit"]][[i]][[j]])
    }
  }
  mean_survfit <- mean_survfit/length(tune_res[["survfit"]])

  mean_conc <- mean_conc[use_index,]
  mean_auc <- mean_auc[use_index,]
  mean_or <- mean_or[use_index,]
  mean_survfit <- mean_survfit[use_index,]

  diff_conc <- as.numeric(test_res[["conc"]][1] - mean_conc[1])
  diff_auc <- test_res[["auc"]][[2]][2] - mean_auc[2]
  diff_or <- test_res[["or"]][2,2] - mean_or[2]
  diff_survfit <- test_res[["survfit"]][nrow(test_res[["survfit"]]),3] - mean_survfit[2]


  ttcomp[k,] <- c(diff_conc, diff_survfit, diff_auc, diff_or)
  k <- k + 1
}


#######################

write_list <- list("authors" = all_author, "method_auc" = method_best_auc, "best_auc" = abs_best_auc, "method_coef" = method_best_coefs, "abs_coef" = abs_best_coefs, "all_coef" = all_poss_coefs, "method_holder" = method_holder, "ttcomp" = ttcomp)
saveRDS(write_list, "all_score_results/meta_results.RDS")
