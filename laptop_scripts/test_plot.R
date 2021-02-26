library(stringr)
library(ggplot2)
library(cowplot)
library(pROC)
library(reshape2)
theme_set(theme_cowplot())

author <- "Shah"
plot_option <- "all"

convert_names <- function(x, the_dict = method_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}

method_names <- read.table("local_info/method_names", stringsAsFactors = F)
disease_namese <- read.table("local_info/disease_names", stringsAsFactors = F, sep = "\t")

all_authors <- unlist(lapply(strsplit(grep("surv", list.files("tune_results/"), value = T), "_", fixed = T), function(x) x[1]))
if(any(all_authors == "Xie")){
  all_authors <- all_authors[-which(all_authors == "Xie")]
}

for(author in all_authors){
  print(author)

res <- readRDS(paste0("test_results/", author, "_res.RDS"))


#CONC ####################################################
conc_df <- data.frame(conc = c(as.numeric(res[["base"]][["conc"]])[1], as.numeric(res[["conc"]])[1]),
                      se = c(as.numeric(res[["base"]][["conc"]])[2], as.numeric(res[["conc"]])[2]),
                      model = c("base", "score"))


conc_plot <- ggplot(conc_df, aes(model, conc)) + geom_point() + 
  geom_errorbar(aes(ymin = conc-1.96*se, ymax = conc+1.96*se), width = 0.2) +
  labs(x = "Model", y = "Concordance") +
  scale_x_discrete(labels = c("Base", "Score"))


#SURVFIT ##################################################
surv_df <- melt(res[["survfit"]][,1:4], id.vars = "time")
se_df <- melt(res[["survfit"]][,c(1, 5:7)], id.vars = "time")
surv_df <- cbind(surv_df, se_df[,-1])
colnames(surv_df) <- c("time", "mean_group", "val_mean", "se_group", "val_se")
surv_df$plot_mean_label <- convert_names( as.character(surv_df$mean_group),
                                          cbind(c("mean_lo", "mean_mid", "mean_hi"), c("Low", "Intermed.", "High")))

surv_plot_1 <- ggplot(surv_df, aes(time/365, val_mean, color = plot_mean_label)) + geom_line() +
  geom_ribbon(aes(ymin = val_mean - val_se, ymax = val_mean + val_se, fill = plot_mean_label),
              alpha = 0.2, color = NA) +
  labs(x = "Years", y = "Cumulative Hazard", color = "Risk\nGroup", fill = "Risk\nGroup") 

surv_df <- surv_df[,-ncol(surv_df)]
base_surv_df <- melt(res[["base"]][["survfit"]][,1:4], id.vars = "time")
se_df <- melt(res[["base"]][["survfit"]][,c(1, 5:7)], id.vars = "time")
base_surv_df <- cbind(base_surv_df, se_df[,-1])
colnames(base_surv_df) <- c("time", "mean_group", "val_mean", "se_group", "val_se")
surv_df$est <- "score"
base_surv_df$est <- "base"

surv_df <- rbind(surv_df, base_surv_df)
surv_plot_2 <- ggplot(surv_df, aes(time/365, val_mean, color = mean_group, linetype = est)) + geom_line() +
  labs(x = "Years", y = "Cumulative Hazard") +
  scale_color_discrete(labels = c("Low", "Intermed.", "High"), name = "Risk\nGroup") +
  scale_linetype_discrete(labels = c("Base", "Score"), name = "Model")

take_time <- min(max(surv_df$time[surv_df$est == "base"]), max(surv_df$time[surv_df$est == "score"]))
sub_surv_df <- surv_df[surv_df$time == take_time,]

surv_plot_3 <- ggplot(sub_surv_df, aes(x = mean_group, y = val_mean, color = est)) + 
  geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.3)) +
  geom_errorbar(aes(ymin = val_mean - val_se, ymax = val_mean + val_se),
                position = position_jitterdodge(jitter.width = 0, dodge.width = 0.3),  width = 0.2) +
  labs(x = "Risk Group", y = "Final Cumulative Hazard") +
  scale_x_discrete(labels = c("Low", "Intermed.", "High")) +
  scale_color_discrete(labels = c("Base", "Score"), name = "Model")


# AUC #########################################################

roc_df=data.frame(tpr = res[["auc"]][[1]]$sensitivities,
                 fpr = 1 - res[["auc"]][[1]]$specificities,
                 model = "score")
base_df=data.frame(tpr = res[["base"]][["auc"]][[1]]$sensitivities,
                  fpr = 1 - res[["base"]][["auc"]][[1]]$specificities,
                  model = "base")
roc_df <- rbind(roc_df, base_df)
auc_plot_1 <- ggplot(roc_df, aes(fpr, tpr, color = model)) + geom_line() +
  geom_abline(intercept=c(0,0),slope=1) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_color_discrete(labels = c("Score", "Base"), name = "Model")


auc_df <- data.frame(rbind(t(matrix(res[["auc"]][[2]])),
                           t(res[["base"]][["auc"]][[2]])))
colnames(auc_df) <- c("lo", "est", "hi")
auc_df$model <- c("score", "base")
p_diff_obj <- roc.test(res[["auc"]][[1]], res[["base"]][["auc"]][[1]])
norm_roc_test_pval <- p_diff_obj$p.value

p_diff <- p_diff_obj$p.value 
if(p_diff < 0.001){
  p_diff <- ">0.001"
} else {
  p_diff <- signif(p_diff, 3)
}


auc_plot_2 <- ggplot(auc_df, aes(model, est)) + geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  labs(x = "Model", y = "AUC", caption = paste0("P-Value = ", p_diff)) +
  scale_x_discrete(labels = c("Base", "Score"))



# ODDS RATIO ###################################################

or_df <- data.frame(res[["or"]])
colnames(or_df) <- c("lo", "val", "hi")
or_df$cut_off <- as.character(c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995))
or_df$model <- "score"
base_df <- data.frame(res[["base"]][["or"]])
colnames(base_df) <- c("lo", "val", "hi")
base_df$cut_off <- as.character(c(0.5, 0.8, 0.9, 0.95, 0.99, 0.995))
base_df$model <- "base"
or_df <- rbind(or_df, base_df)

or_plot <- ggplot(or_df, aes(cut_off, val, color = model)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.3)) +
  geom_errorbar(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.3), 
                aes(ymin = lo, ymax = hi), width = 0.2) +
  labs(x = "Cut Off", y = "Odds Ratio") +
  scale_color_discrete(labels = c("Base", "Score"), name = "Model")



# PR ###################################################

pr_df=data.frame(recall = res[["pr"]]$curve[,1],
                  precision = res[["pr"]]$curve[,2],
                  model = "score")
base_df=data.frame(recall = res[["base"]][["pr"]]$curve[,1],
                   precision = res[["base"]][["pr"]]$curve[,2],
                   model = "base")
pr_df <- rbind(pr_df, base_df)
the_cap <- paste0("AUC: Base = ", signif(res[["base"]][["pr"]]$auc.integral, 3), ", Score = ", 
                  signif(res[["pr"]]$auc.integral, 3))
pr_plot <- ggplot(pr_df, aes(recall, precision, color = model)) + geom_line() +
  labs(x = "Recall", y = "Precision", caption = the_cap) +
  scale_color_discrete(labels = c("Score", "Base"), name = "Model")
pr_auc <- data.frame(pr_auc = c(res[["base"]][["pr"]]$auc.integral, res[["pr"]]$auc.integral),
                     model = c("base", "score"))

## OTHER SCORES #################################################

all_other_plots <- list()
ll <- 1
if(!is.null(res[["other"]][[1]])){
  author_defs <- read.table("test_results/author_scores", stringsAsFactors = F, header = T)
  author_defs <- author_defs[author_defs$Author == author,]
  author_defs <- author_defs[order(as.numeric(substr(author_defs$PGS_Catalog_name, 4, nchar(author_defs$PGS_Catalog_name)))),]
  
  #all_other_scores <- readRDS("test_results/all_score.1.RDS")
  #author_defs <- author_defs[author_defs$PGS_Catalog_name %in% colnames(all_other_scores),]
  
  all_other_auc <- list()
  all_other_roc_test <- list()
  if(class(res[["other"]][["auc"]][[1]]) == "roc"){
    iter_len <- 1
  } else {
    iter_len <- length(res[["other"]][["auc"]][[1]])
  }
  for(kk in 1:iter_len){
    all_other_roc_test[[kk]] <- roc.test(res[["other"]][["auc"]][[1]][[1]], res[["auc"]][[1]])$p.value
    other_auc <- as.data.frame(rbind(res[["auc"]][[2]],
                       res[["other"]][["auc"]][[2]][[kk]]))
    colnames(other_auc) <- c("lo", "auc", "hi")
    other_auc$type <- c("Internal", paste(author_defs$Author.1[kk], author_defs$Year[kk], sep = "-"))
    all_other_auc[[kk]] <- other_auc
    
    all_other_plots[[ll]] <- ggplot(other_auc, aes(auc, type)) + geom_point() +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
      labs(x = "AUC", y = "Type")
  }
} else{
  all_other_auc <- NULL
  all_other_roc_test <- NULL
}

## EXTRA COVARS #####################################################

if(!is.null(res[["extra"]][[1]])){
  extra_all_roc_test <- roc.test(res[["extra"]][["base_auc"]][[1]], res[["extra"]][["auc"]][[1]])$p.value
  #res[["extra"]][["auc"]][[1]]
  #res[["extra"]][["base_auc"]][[1]]
  #res[["auc"]][[1]]
  #res[["base"]][["auc"]][[1]]
  
  extra_auc <- as.data.frame(rbind(res[["auc"]][[2]], 
    res[["base"]][["auc"]][[2]],
    res[["extra"]][["auc"]][[2]],
    res[["extra"]][["base_auc"]][[2]]))
  colnames(extra_auc) <- c("lo", "auc", "hi")
  extra_auc$type <- c("Basic", "Basic", "Extra", "Extra")
  extra_auc$Score <- c("Yes", "No", "Yes", "No")
  
  extra_plot <- ggplot(extra_auc, aes(auc, type, color = Score)) + geom_point() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.1) +
    labs(x = "AUC", y = "Covariate Type", color = "Include\nScore") 

} else {
  extra_plot <- NULL
  extra_auc <- NULL
  extra_all_roc_test <- NULL
}


###################### to save ###########################

to_save <- list(conc_df, surv_df, sub_surv_df, roc_df, auc_df, pr_auc, or_df, all_other_auc, extra_auc)
saveRDS(to_save, paste0("test_results/", author, ".best_data.RDS"))
roc_save <- list(norm_roc_test_pval, all_other_roc_test, extra_all_roc_test)
saveRDS(roc_save, paste0("test_results/", author, ".roc_data.RDS"))

# PLOT ###################################################

if(plot_option == "all"){
  ggsave(paste0("test_plots/", tolower(author), ".conc.test.png"),
         conc_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".survfit1.test.png"),
         surv_plot_1 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".survfit2.test.png"),
         surv_plot_2 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".survfit3.test.png"),
         surv_plot_3 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".auc1.test.png"),
         auc_plot_1 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".auc2.test.png"),
         auc_plot_2 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".or.test.png"),
         or_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("test_plots/", tolower(author), ".pr.test.png"),
         pr_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  if(length(all_other_plots) > 0){
    for(kk in 1:length(all_other_plots)){
      ggsave(paste0("test_plots/", tolower(author), ".other.test.", kk, ".png"),
         all_other_plots[[kk]] + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    }
  }
  
  if(!is.null(extra_plot)){
    ggsave(paste0("test_plots/", tolower(author), ".extra.test.png"),
         extra_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  }
}

}

