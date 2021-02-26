library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

author <- "Kottgen"

res <- readRDS(paste0("tune_results/", author, "_res.RDS"))

split_score_names <- str_split(res[["score_names"]], fixed("."), simplify = T)
simple_name <- paste0(split_score_names[,3], "-", split_score_names[,2])
method_name <- split_score_names[,3]


#do_plot_series <- function(subset_inds, draw_plot_option = "none", save_option = FALSE, best_stat = "auc"){
  
  if(length(unique(method_name[subset_inds])) == 1){
    plot_ext <- unique(method_name[subset_inds])
  } else {
    plot_ext <- "best"
  }
  
  ##############################################################################
  #Average all of the results
  mean_conc <- Reduce("+", res[["conc"]])/length(res[["conc"]])
  
  mean_survfit <- matrix(0, nrow = length(res[["survfit"]][[1]]), ncol = 6)
  for(i in 1:length(res[["survfit"]])){
    for(j in 1:length(res[["score_names"]])){
      mean_survfit[j,] <- mean_survfit[j,] + as.numeric(tail(res[["survfit"]][[i]][[j]],1)[-1])
    }
  }
  mean_survfit <- mean_survfit/length(res[["survfit"]])
  
  mean_auc <- Reduce("+", res[["auc"]])/length(res[["auc"]])
  
  mean_or <- Reduce("+", res[["or"]])/length(res[["or"]])
  
  
  ###########################################################################
  #Average for the base items
  base_conc <- rep(0, 3)
  for(i in 1:length(res[["base"]][[1]])){
    base_conc <- base_conc + as.numeric(c(res[["base"]][[1]][[i]]$concordance,
                                          res[["base"]][[1]][[i]]$std.err,
                                          res[["base"]][[1]][[i]]$std.err * 1.96))
  }
  base_conc <- base_conc/length(res[["base"]][[1]])
  base_conc_mat <- data.frame("base", "base", t(base_conc), stringsAsFactors = F)
  colnames(base_conc_mat) <- c("score_name", "method", "conc", "se", "ci")
  
  base_survfit <- matrix(0, nrow = 1, ncol = 6)
  for(i in 1:length(res[["base"]][[2]])){
    base_survfit <- base_survfit + as.numeric(tail(res[["base"]][[2]][[i]],1)[-1])
  }
  base_survfit <- base_survfit/9
  
  base_auc <- 0
  for(i in 1:length(res[["base"]][[3]])){
    base_auc <- res[["base"]][[3]][[i]] + base_auc
  }
  base_auc <- base_auc/length(res[["base"]][[3]])
  
  base_or <- 0
  for(i in 1:length(res[["base"]][[4]])){
    base_or <- res[["base"]][[4]][[i]] + base_or
  }
  base_or <- base_or/length(res[["base"]][[4]])


