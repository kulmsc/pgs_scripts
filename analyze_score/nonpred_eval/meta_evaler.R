
#When all of the tuning is done want to go though and check that 
#If method performance is dependent on any factors (sample size, SNPs, etc.)
  #Do for one AUC for each author (the best)
  #Do this once for each method (the best method for each author)

#Want to see if there is a correlation between method parameters and the meta stats
  #so collect all results (for all authors) then run a model

#Also want to compare training to testing results

#Will focus on AUC

list_all_auc <- list()
list_all_best_names <- list()
list_all_best_methods <- list()
list_ukb_cc <- list()
list_ukb_snps <- list()

i <- 1
all_author <- unlist(lapply(strsplit(list.files("../tune_score/tune_results/", "res"), "_"), function(x) x[1]))
for(author in all_author){
  all_res <- readRDS(paste0("../tune_score/tune_results/", author, "_res.RDS"))
  score_names <- all_res[["score_names"]]

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
  ukb_snps <- rep(0, length(score_names))
  for(sname in score_names){
    new_name <- paste0(strsplit(sname, ".", fixed = T)[[1]][3], ".", strsplit(sname, ".", fixed = T)[[1]][2], ".ss")
    get_files <- list.files(paste0("../../mod_sets/", author, "/"), pattern = new_name)
    temp_length <- rep(0, length(get_files))
    for(j in 1:length(get_files)){
      temp_length[j] <- as.numeric(system(paste0("cat ../../mod_sets/", author, "/", get_files[j], " | wc -l"), intern = TRUE))
    }
    ukb_snps[which(score_names == sname)] <- sum(temp_length)
  }

  names(ukb_snps) <- score_names
  list_ukb_snps[[i]] <- ukb_snps

  list_all_auc[[i]] <- mean_auc
  i <- i + 1
}

#add meta info to the all_auc
all_auc <- do.call("rbind", list_all_auc)
all_auc <- cbind(all_auc, do.call("rbind", list_ukb_cc), unlist(list_ukb_snps))
colnames(all_auc) <- c("auc_lo", "auc", "auc_hi", "name", "method", "author", "uk_control", "uk_case", "uk_snps")

meta_info <- read.table("../../raw_ss/meta_stats", stringsAsFactors=F, sep=",", header=T)
add_on <- data.frame(matrix(0, nrow = nrow(all_auc), ncol = 10))
colnames(add_on) <- colnames(meta_info)[-c(1,2)]
for(author in unique(all_auc$author)){
  add_on[all_auc$author == author,] <- meta_info[meta_info$author == author,-c(1,2)]
}
all_auc <- cbind(all_auc, add_on)

#get the best names (absolute best, and best name for each method for each author)
abs_best_names <- unlist(lapply(list_all_best_names, function(x) x[3,2]))
method_best_names <- unlist(list_all_best_methods)

#subset the all_auc
abs_best_auc <- all_auc[all_auc$name %in% abs_best_names,]
method_best_auc <- all_auc[all_auc$name %in% method_best_names,]


############################################################
################ META INFO REGRESSIONS #####################
############################################################

meta_to_regress <- c("sampe_size", "cases", "cases/controls", "snps", "h2", "uk_snps", "uk_case/uk_control")
abs_best_coefs <- matrix(0, nrow = length(meta_to_regress), ncol = 4)
method_best_coefs <- rep(list(matrix(0, nrow = length(meta_to_regress), ncol = 4)), length(unique(method_best_auc$method)))

for(i in 1:length(meta_to_regress)){
  best_mod <- lm(as.formula(paste0("auc ~ ", meta_to_regress[i])), data = abs_best_auc)
  abs_best_coefs[i,] <- summary(best_mod)$coefficients[2,]

  for(m in unique(method_best_auc$method)){
     method_mod <- lm(as.formula(paste0("auc ~ ", meta_to_regress[i])), data = method_best_auc[method_best_auc$method == m,])
     method_best_coefs[[which(unique(method_best_auc$method) == m)]][i,] <- summary(method_mod)$coefficients[2,]
  }
}


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
      mean_survfit[j,] <- mean_survfit[j,] + as.numeric(tail(tune_res[["survfit"]][[i]][[j]],1)[-1])
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

write_list <- list("authors" = all_author, "method_auc" = method_best_auc, "best_auc" = abs_best_auc, "method_coef" = method_best_coefs, "abs_coef" = abs_best_coefs, "ttcomp" = ttcomp)
saveRDS(write_list, "all_score_results/meta_results.RDS")
