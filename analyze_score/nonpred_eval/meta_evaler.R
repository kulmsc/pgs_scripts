
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
  list_ukb_cc <- as.numeric(table(rowSums(ukb_pheno[,1:4]) > 0))

  #get the ukb snp lengths
  ukb_snps <- rep(0, length(score_names))
  for(sname in score_names){
    get_files <- list.files("../../mod_sets/IMSGC/", pattern = "imsgc.1.clump")
    temp_length <- rep(0, length(get_files))
    for(j in 1:length(get_files)){
      temp_length[j] <- as.numeric(system(paste0("cat ../../mod_sets/", author, "/", get_files[j], " | wc -l")))
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

meta_to_regress <- c("sampe_size", "cases", "cases/controls", "snps", "h2")
abs_best_model_list <- list()
for(i in 1:length(meta_to_regress)){
  lm(V2 + 
}
