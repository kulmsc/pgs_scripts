library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

author <- "Christophersen"
plot_option <- "all"

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


surv_plot_1 <- ggplot(surv_df, aes(time/365, val_mean, color = mean_group)) + geom_line() +
  geom_ribbon(aes(ymin = val_mean - val_se, ymax = val_mean + val_se, fill = mean_group),
              alpha = 0.2, color = NA) +
  labs(x = "Years", y = "Cumulative Hazard")

base_surv_df <- melt(res[["base"]][["survfit"]][,1:4], id.vars = "time")
se_df <- melt(res[["base"]][["survfit"]][,c(1, 5:7)], id.vars = "time")
base_surv_df <- cbind(base_surv_df, se_df[,-1])
colnames(base_surv_df) <- c("time", "mean_group", "val_mean", "se_group", "val_se")
surv_df$est <- "score"
base_surv_df$est <- "base"

surv_df <- rbind(surv_df, base_surv_df)
surv_plot_2 <- ggplot(surv_df, aes(time/365, val_mean, color = mean_group, linetype = est)) + geom_line() +
  labs(x = "Years", y = "Cumulative Hazard") +
  scale_color_discrete(labels = c("Low", "Mid", "High"), name = "Risk") +
  scale_linetype_discrete(labels = c("Base", "Score"), name = "Model")

take_time <- min(max(surv_df$time[surv_df$est == "base"]), max(surv_df$time[surv_df$est == "score"]))
sub_surv_df <- surv_df[surv_df$time == take_time,]

surv_plot_3 <- ggplot(sub_surv_df, aes(x = mean_group, y = val_mean, color = est)) + geom_point() +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = val_mean - val_se, ymax = val_mean + val_se),
                position = position_dodge(width = 0.3),  width = 0.2) +
  labs(x = "Risk", y = "Cumulative Hazard") +
  scale_x_discrete(labels = c("Low", "Mid", "High")) +
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
p_diff <- roc.test(res[["auc"]][[1]], res[["base"]][["auc"]][[1]])$p
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

or_plot <- ggplot(or_df, aes(cut_off, val, color = model)) + geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
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
  









# PLOT ###################################################

if(plot_option == "all"){
  ggsave(paste0("output_plots/", tolower(author), ".conc.test.png"),
         conc_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".survfit1.test.png"),
         surv_plot_1 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".survfit2.test.png"),
         surv_plot_2 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".survfit3.test.png"),
         surv_plot_3 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".auc1.test.png"),
         auc_plot_1 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".auc2.test.png"),
         auc_plot_2 + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".or.test.png"),
         or_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
  
  ggsave(paste0("output_plots/", tolower(author), ".pr.test.png"),
         pr_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
}


