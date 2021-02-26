library(stringr)
library(ggplot2)
library(cowplot)
library(viridis)
library(stringr)
theme_set(theme_cowplot())



all_best_names <- list.files("test_results/", "best")
best_author <-  str_split(all_best_names, fixed("."), simplify = T)[,1]

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
author_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = '\t')

convert_names <- function(x, the_dict = author_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}


all_auc <- list()
all_conc <- list()
all_or <- list()
all_km <- list()
all_subkm <- list()
all_roc <- list()
all_pr <- list()
all_other_auc <- list()
all_extra_auc <- list()

for(i in 1:length(best_author)){
  curr_best <- readRDS(paste0("test_results/", all_best_names[i]))
  
  all_conc[[i]] <- curr_best[[1]]
  all_conc[[i]] <- rbind(all_conc[[i]], data.frame(all_conc[[i]][2,1:2] - all_conc[[i]][1,1:2], model = "diff"))
  all_conc[[i]]$author <- best_author[i]
  
  all_km[[i]] <- curr_best[[2]]
  all_km[[i]]$author <- best_author[i]
  all_km[[i]]$diff_val <- 0
  all_km[[i]]$diff_val[all_km[[i]]$mean_group == "mean_lo"] <- all_km[[i]]$val_mean[all_km[[i]]$mean_group == "mean_lo"] -
    all_km[[i]]$val_mean[all_km[[i]]$mean_group == "mean_mid"]
  all_km[[i]]$diff_val[all_km[[i]]$mean_group == "mean_hi"] <- all_km[[i]]$val_mean[all_km[[i]]$mean_group == "mean_hi"] -
    all_km[[i]]$val_mean[all_km[[i]]$mean_group == "mean_mid"]
  
  all_subkm[[i]] <- curr_best[[3]]
  all_subkm[[i]]$author <- best_author[i]
  
  all_roc[[i]] <- curr_best[[4]]
  all_roc[[i]]$author <- best_author[i]
  
  all_auc[[i]] <- curr_best[[5]]
  all_auc[[i]] <- rbind(all_auc[[i]], data.frame(all_auc[[i]][1,1:3] - all_auc[[i]][2,1:3], model = "diff"))
  all_auc[[i]]$author <- best_author[i]
  
  all_pr[[i]] <- curr_best[[6]]
  all_pr[[i]] <- rbind(all_pr[[i]], data.frame(pr_auc = all_pr[[i]][2,1] - all_pr[[i]][1,1], model = "diff"))
  all_pr[[i]]$author <- best_author[i]
  
  all_or[[i]] <- curr_best[[7]]
  curr_best[[7]]$val[curr_best[[7]]$model == "score"] / curr_best[[7]]$val[curr_best[[7]]$model == "score"][1]
  all_or[[i]] <- rbind(all_or[[i]], data.frame(all_or[[i]][all_or[[i]]$model == "score",1:3] - all_or[[i]][all_or[[i]]$model == "base",1:3],
             cut_off = unique(all_or[[i]]$cut_off), model = "diff"))
  all_or[[i]]$author <- best_author[i]
  
  all_other_auc[[i]] <- curr_best[[8]]
  if(!is.null(curr_best[[8]])){
    all_other_auc[[i]]$author <- best_author[i]
  }
  
  all_extra_auc[[i]] <- curr_best[[9]]
  if(!is.null(curr_best[[9]])){
    all_extra_auc[[i]]$author <- best_author[i]
  }
}



################################################################################################
#                                 AUC                                  ####
#################################################################################################

keep_all_auc <- all_auc
all_auc <- do.call("rbind", all_auc)
all_auc$disease <- convert_names(all_auc$author)

all_abs_auc <- all_auc[all_auc$model %in% c("base", "score"),]
ranked_disease <- all_abs_auc$disease[all_abs_auc$model == "score"][order(all_abs_auc$est[all_abs_auc$model == "score"])]
all_abs_auc$disease <- factor(all_abs_auc$disease, levels = ranked_disease)
all_abs_auc$model <- str_to_title(all_abs_auc$model)

the_plot <- ggplot(all_abs_auc, aes(est, disease, color = model)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, 
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
  geom_hline(yintercept = 1:length(unique(all_abs_auc$disease)) + 0.5, color = "grey80") +
  labs(x = "AUC", color = "Model", y = "")
plot(the_plot)
ggsave(paste0("test_plots/meta/auc.val.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 6.5, height = 5)

table_col_1 <- all_abs_auc[all_abs_auc$model == "Score",]

all_diff_auc <- all_auc[all_auc$model == "diff",]
all_diff_auc$disease <- factor(all_diff_auc$disease, levels = all_diff_auc$disease[order(all_diff_auc$est)])

the_plot <- ggplot(all_diff_auc, aes(est, disease)) +
  geom_point() + labs(x = "AUC Improvement", y = "") +
  geom_hline(yintercept = 1:length(unique(all_diff_auc$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("test_plots/meta/auc.diff.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 5)

table_col_2 <- all_diff_auc

#stats
min(all_auc$est[all_auc$model == "score"])
max(all_auc$est[all_auc$model == "score"])
all_auc$disease[all_auc$est == max(all_auc$est[all_auc$model == "score"])]
all_auc$disease[all_auc$est == min(all_auc$est[all_auc$model == "score"])]

min(all_auc$est[all_auc$model == "diff"])
max(all_auc$est[all_auc$model == "diff"])
all_auc$disease[all_auc$est == max(all_auc$est[all_auc$model == "diff"])]
all_auc$disease[all_auc$est == min(all_auc$est[all_auc$model == "diff"])]

cor(all_auc$est[all_auc$model == "score"], all_auc$est[all_auc$model == "diff"], method = "spearman")


################################################################################################
#                                 ODDS RATIO                                ####
#################################################################################################

all_or <- do.call("rbind", all_or)
all_or$disease <- convert_names(all_or$author)
all_or$val[all_or$val > 100] <- 100
all_or$hi[all_or$hi > 100] <- 100
all_or$lo[all_or$lo > 100] <- 100

#One plot for each cut off
for(curr_cu in unique(all_or$cut_off)){
    sub_or <- all_or[all_or$cut_off == curr_cu & all_or$model %in% c("score", "base"),]
    ranked_disease <- sub_or$disease[sub_or$model == "score"][order(sub_or$val[sub_or$model == "score"])]
    sub_or$disease <- factor(sub_or$disease, levels = ranked_disease)
    sub_or$model <- str_to_title(sub_or$model)
    
    the_plot <- ggplot(sub_or, aes(x = val, y = disease, color = model)) + 
      geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, 
                     position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
      labs(x = "Odds Ratio", y = "", color = "Model", caption = paste0("Cutoff = ", curr_cu)) +
      geom_hline(yintercept = 1:length(unique(all_diff_auc$disease)) + 0.5, color = "grey80")
    plot(the_plot)
    ggsave(paste0("test_plots/meta/or.", curr_cu, ".png"),
           the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 6)
}

#Absolute scores
sub_or <- all_or[all_or$model == "score",]
ranked_disease <- sub_or$disease[sub_or$cut_off == "0.9"][order(sub_or$val[sub_or$cut_off == "0.9"])]
sub_or$disease <- factor(sub_or$disease, levels = ranked_disease)

top_names <- as.character(unique(sub_or$disease[sub_or$val > 25]))
bottom_names <- as.character(unique(sub_or$disease))[!as.character(unique(sub_or$disease)) %in% top_names]


top_plot <- ggplot(sub_or[sub_or$disease %in% top_names,], aes(x = val, y = disease, color = cut_off)) + geom_point() +
  labs(x = "", y = "", color = "Cut Off") +
  geom_hline(yintercept = 1:length(top_names) + 0.5, color = "grey80")
bottom_plot <- ggplot(sub_or[sub_or$disease %in% bottom_names,], aes(x = val, y = disease, color = cut_off)) + geom_point() +
  labs(x = "Odds Ratio", y = "", color = "Cut Off") +
  geom_hline(yintercept = 1:length(bottom_names) + 0.5, color = "grey80")
the_plot <- plot_grid(top_plot, bottom_plot, ncol = 1, rel_heights = c(1, 2.7))
plot(the_plot)

ggsave(paste0("test_plots/meta/or.val.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 6)
ggsave(paste0("test_plots/meta/or.valtop.png"),
       top_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 3.5)
ggsave(paste0("test_plots/meta/or.valbottom.png"),
       bottom_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 5)

table_col_3 <- sub_or[sub_or$cut_off == 0.95,]

#Difference
sub_or <- all_or[all_or$model == "diff",]
ranked_disease <- sub_or$disease[sub_or$cut_off == "0.9"][order(sub_or$val[sub_or$cut_off == "0.9"])]
sub_or$disease <- factor(sub_or$disease, levels = ranked_disease)

top_names <- as.character(unique(sub_or$disease[sub_or$val > 25]))
bottom_names <- as.character(unique(sub_or$disease))[!as.character(unique(sub_or$disease)) %in% top_names]


top_plot <- ggplot(sub_or[sub_or$disease %in% top_names,], aes(x = val, y = disease, color = cut_off)) + geom_point() +
  labs(x = "", y = "", color = "Cut Off") +
  geom_hline(yintercept = 1:length(top_names) + 0.5, color = "grey80")
bottom_plot <- ggplot(sub_or[sub_or$disease %in% bottom_names,], aes(x = val, y = disease, color = cut_off)) + geom_point() +
  labs(x = "Odds Ratio", y = "", color = "Cut Off") +
  geom_hline(yintercept = 1:length(bottom_names) + 0.5, color = "grey80")
the_plot <- plot_grid(top_plot, bottom_plot, ncol = 1, rel_heights = c(1, 2))
plot(the_plot)

ggsave(paste0("test_plots/meta/or.diff.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 6)
ggsave(paste0("test_plots/meta/or.difftop.png"),
       top_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 3.5)
ggsave(paste0("test_plots/meta/or.diffbottom.png"),
       bottom_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 5)


table_col_4 <- sub_or[sub_or$cut_off == 0.95,]

total_table <- data.frame(author = table_col_1$disease[order(table_col_1$author)],
                          auc_total = table_col_1$est[order(table_col_1$author)],
                          auc_diff = table_col_2$est[order(table_col_2$author)],
                          or_total = table_col_3$val[order(table_col_3$author)],
                          or_diff = table_col_4$val[order(table_col_4$author)])
write.table(total_table, "test_results/for_table.txt", row.names = F, col.names = T, quote = F, sep = "\t")


#stats
or_95_score <- all_or[all_or$cut_off == 0.95 & all_or$model == "score",]
or_95_base <- all_or[all_or$cut_off == 0.95 & all_or$model == "base",]

use_or <- all_or[all_or$cut_off == 0.95,]

min(use_or$val[use_or$model == "score"])
max(use_or$val[use_or$model == "score"])
use_or$disease[use_or$val == max(use_or$val[use_or$model == "score"])]
use_or$disease[use_or$val == min(use_or$val[use_or$model == "score"])]

min(use_or$val[use_or$model == "diff"])
max(use_or$val[use_or$model == "diff"])
use_or$disease[use_or$val == max(use_or$val[use_or$model == "diff"])]
use_or$disease[use_or$val == min(use_or$val[use_or$model == "diff"])]

exit()
cor(use_or$val[use_or$model == "score"], use_or$val[all_auc$model == "diff"], method = "spearman")

all_or[all_or$disease == "Prostate Cancer" & all_or$cut_off == 0.995,]
all_or[all_or$disease == "A. Fib." & all_or$cut_off == 0.995,]

################################################################################################
#                                 CONCORDANCE                               ####
#################################################################################################

all_conc <- do.call("rbind", all_conc)
all_conc$disease <- convert_names(all_conc$author)


all_abs_conc <- all_conc[all_conc$model %in% c("base", "score"),]
ranked_disease <- all_abs_conc$disease[all_abs_conc$model == "score"][order(all_abs_conc$conc[all_abs_conc$model == "score"])]
all_abs_conc$disease <- factor(all_abs_conc$disease, levels = ranked_disease)
all_conc$model <- str_to_title(all_conc$model)

the_plot <- ggplot(all_abs_conc, aes(conc, disease, color = model)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
  geom_errorbarh(aes(xmin = conc - se, xmax = conc + se), height = 0, 
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
  geom_hline(yintercept = 1:length(unique(all_abs_conc$disease)) + 0.5, color = "grey80") +
  labs(x = "Concordance", color = "Model", y = "")
plot(the_plot)
ggsave(paste0("test_plots/meta/conc.val.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 6)

all_diff_conc <- all_conc[all_conc$model == "diff",]
all_diff_conc$disease <- factor(all_diff_conc$disease, levels = all_diff_conc$disease[order(all_diff_conc$conc)])
all_conc$model <- str_to_title(all_conc$model)

the_plot <- ggplot(all_diff_conc, aes(conc, disease)) +
  geom_point() + labs(x = "Concordance Improvement", y = "") +
  geom_hline(yintercept = 1:length(unique(all_diff_conc$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("test_plots/meta/conc.diff.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 6)


################################################################################################
#                                 KM                            ####
#################################################################################################

#More advanced comparisons that require substracting from base and dividing between groups is likely too complex to be worth it

all_km <- do.call("rbind", all_km)
all_km$disease <- convert_names(all_km$author)

diff_subkm <- list()
for(i in 1:length(all_subkm)){
  diff_subkm[[i]] <- all_subkm[[i]][all_subkm[[i]]$est == "score", ]
  diff_subkm[[i]]$val_mean <- all_subkm[[i]]$val_mean[all_subkm[[i]]$est == "score"] -
    all_subkm[[i]]$val_mean[all_subkm[[i]]$est == "base"]
  diff_subkm[[i]]$est <- "diff"
}

all_subkm <- do.call("rbind", all_subkm)
all_subkm$disease <- convert_names(all_subkm$author)

diff_subkm <- do.call("rbind", diff_subkm)
diff_subkm$disease <- convert_names(diff_subkm$author)

dubsub_km <- all_subkm[all_subkm$mean_group == "mean_hi" & all_subkm$est == "score",]
order_author <- dubsub_km$author[order(dubsub_km$val_mean)]
start_inds <- c(1, 5, 9, 13, 17, 21)
stop_inds <- c(4, 8, 12, 16, 20, 24)

all_sub_plots <- list()
for(i in 1:6){
  plot_subkm <- all_km[all_km$author %in% order_author[start_inds[i]:stop_inds[i]] & all_km$est == "score",]
  all_sub_plots[[i]] <- ggplot(plot_subkm, aes(time/365, diff_val, color = disease,
                                               group = interaction(mean_group, disease))) + 
    geom_line(size = 1.5) + labs(x = "Years", y = "Cumulative Hazard", color = "")
}
the_plot <- plot_grid(plotlist = all_sub_plots, nrow = 2)
plot(the_plot)
ggsave(paste0("test_plots/meta/km.all.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 15, height = 9.5)


plot_subkm <- all_subkm[all_subkm$est == "score",]
disease_order <- plot_subkm$disease[plot_subkm$mean_group == "mean_mid"][order(plot_subkm$val_mean[plot_subkm$mean_group == "mean_mid"])]
plot_subkm$disease <- factor(plot_subkm$disease, levels = disease_order)
the_plot <- ggplot(plot_subkm, aes(val_mean, disease, color = mean_group)) + geom_point() +
  geom_hline(yintercept = 1:length(unique(plot_subkm$disease)) + 0.5, color = "grey80") +
  labs(x = "Final Cumulative Incidence", y = "", color = "Risk Group") +
  scale_color_discrete(labels = c("Low", "Intermed.", "High"))
plot(the_plot)
ggsave(paste0("test_plots/meta/kmfinal.val.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 6, height = 6)

disease_order <- diff_subkm$disease[diff_subkm$mean_group == "mean_mid"][order(diff_subkm$val_mean[diff_subkm$mean_group == "mean_mid"])]
diff_subkm$disease <- factor(diff_subkm$disease, levels = rev(disease_order))
the_plot <- ggplot(diff_subkm, aes(val_mean, disease, color = mean_group)) + geom_point() +
  geom_hline(yintercept = 1:length(unique(diff_subkm$disease)) + 0.5, color = "grey80") +
  geom_vline(xintercept = 0) +
  labs(x = "Difference in Final Cumulative Hazard", y = "", color = "Risk Group") +
  scale_color_discrete(labels = c("Low", "Intermed.", "High"))
plot(the_plot)
ggsave(paste0("test_plots/meta/kmfinal.diff.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 6, height = 6)


##############################################################################################
for_corr <- data.frame("auc" = all_auc$est[all_auc$model == "score"],
                       "aucimp" = all_auc$est[all_auc$model == "diff"],
                       "or" = all_or$val[all_or$model == "score" & all_or$cut_off == 0.95],
                       "orimp" = all_or$val[all_or$model == "diff" & all_or$cut_off == 0.95],
                       "conc" = all_conc$conc[all_conc$model == "Score"],
                       "concimp" = all_conc$conc[all_conc$model == "Diff"],
                       "km" = plot_subkm$val_mean[plot_subkm$mean_group == "mean_hi" & plot_subkm$est == "score"]/
                              plot_subkm$val_mean[plot_subkm$mean_group == "mean_lo" & plot_subkm$est == "score"],
                       "diffkm"  =  diff_subkm$val_mean[diff_subkm$mean_group == "mean_hi" & diff_subkm$est == "diff"]/
                              diff_subkm$val_mean[diff_subkm$mean_group == "mean_lo" & diff_subkm$est == "diff"])
cor_df <- signif(cor(for_corr, method = "spearman"), 3)
write.table(cor_df, "supp_tables/cor_test.txt", row.names = T, col.names = T, quote = F, sep = "\t")

write_subkm <- data.frame("disease" = plot_subkm$disease[plot_subkm$mean_group == "mean_lo"],
                          "lo" = signif(plot_subkm$val_mean[plot_subkm$mean_group == "mean_lo"], 3),
                          "mid" = signif(plot_subkm$val_mean[plot_subkm$mean_group == "mean_mid"], 3),
                          "hi" = signif(plot_subkm$val_mean[plot_subkm$mean_group == "mean_hi"], 3))
write.table(write_subkm, "supp_tables/subkm_test.txt", row.names = F, col.names = T, quote = F, sep = "\t")

write_diffkm <- data.frame("disease" = diff_subkm$disease[diff_subkm$mean_group == "mean_lo"],
                          "lo" = signif(diff_subkm$val_mean[diff_subkm$mean_group == "mean_lo"], 3),
                          "mid" = signif(diff_subkm$val_mean[diff_subkm$mean_group == "mean_mid"], 3),
                          "hi" = signif(diff_subkm$val_mean[diff_subkm$mean_group == "mean_hi"], 3))
write.table(write_diffkm, "supp_tables/diffkm_test.txt", row.names = F, col.names = T, quote = F, sep = "\t")


                          

################################################################################################
#                                 PR                            ####
#################################################################################################

all_pr <- do.call("rbind", all_pr)
all_pr$disease <- convert_names(all_pr$author)

all_abs_auc <- all_pr[all_pr$model %in% c("base", "score"),]
ranked_disease <- all_abs_auc$disease[all_abs_auc$model == "score"][order(all_abs_auc$pr_auc[all_abs_auc$model == "score"])]
all_abs_auc$disease <- factor(all_abs_auc$disease, levels = ranked_disease)
all_abs_auc$model <- str_to_title(all_abs_auc$model)

the_plot <- ggplot(all_abs_auc, aes(pr_auc, disease, color = model)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.3)) +
  geom_hline(yintercept = 1:length(unique(all_abs_auc$disease)) + 0.5, color = "grey80") +
  labs(x = "PR AUC", color = "Model", y = "")
plot(the_plot)
ggsave(paste0("test_plots/meta/pr.val.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5.5, height = 6)


all_diff_auc <- all_pr[all_pr$model == "diff",]
all_diff_auc$disease <- factor(all_diff_auc$disease, levels = all_diff_auc$disease[order(all_diff_auc$pr_auc)])
all_diff_auc$model <- str_to_title(all_diff_auc$model)

the_plot <- ggplot(all_diff_auc, aes(pr_auc, disease)) +
  geom_point() + labs(x = "PR AUC Improvement", y = "") +
  geom_hline(yintercept = 1:length(unique(all_diff_auc$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("test_plots/meta/pr.diff.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 6, height = 6)


################################################################################################
#                                 all other score                            ####
#################################################################################################

j <- 1
better_all_other <- list()
for(i in 1:length(all_other_auc)){
  if(!is.null(all_other_auc[[i]])){
    temp_list <- list()
    for(k in 1:(length(all_other_auc[[i]])-1)){
      if(k == 1){
        temp_list[[k]] <- all_other_auc[[i]][[k]]
      } else {
        temp_list[[k]] <- all_other_auc[[i]][[k]][2,]
      }
    }
    better_all_other[[j]] <- do.call("rbind", temp_list)
    better_all_other[[j]]$author <- all_other_auc[[i]][[length(all_other_auc[[i]])]]
    better_all_other[[j]]$source <- c("Internal", rep("External", nrow(better_all_other[[j]])-1))
    print(i)
    j <- j + 1
  }
}

quick_fun <- function(x){
  x$auc[x$source == "Internal"] - max(x$auc[x$source == "External"])
}

mean(unlist(lapply(better_all_other, quick_fun)))
sd(unlist(lapply(better_all_other, quick_fun)))/sqrt(length(better_all_other))


better_all_other <- do.call("rbind", better_all_other)
better_all_other$disease <- convert_names(better_all_other$author)
auth_order <- better_all_other$disease[better_all_other$source == "Internal"][
  order(better_all_other$auc[better_all_other$source == "Internal"])]
better_all_other$disease <- factor(better_all_other$disease, levels = auth_order)
the_plot <- ggplot(better_all_other, aes(auc, disease, color = source)) + geom_point() +
  labs(x = "AUC", y = "", color = "Score\nSource") +
  geom_hline(yintercept = 1:length(unique(better_all_other$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("test_plots/meta/other.score.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 6.5, height = 5)

write_other <- data.frame("disease" = better_all_other$disease,
                          "AUC" = signif(better_all_other$auc, 3),
                          "type" = better_all_other$type)
write.table(write_other, "supp_tables/other_score_test.txt", row.names = F, col.names = T, quote = F, sep = "\t")

################################################################################################
#                                 all other covars                         ####
#################################################################################################

better_all_extra <- list()
for(i in 1:length(all_extra_auc)){
  if(!is.null(all_extra_auc[[i]])){
    better_all_extra[[i]] <- all_extra_auc[[i]][all_extra_auc[[i]]$type == "Extra",][,-4]
    colnames(better_all_extra[[i]]) <- c("lo", "est", "hi", "model", "author")
    better_all_extra[[i]]$model <- c("score", "base")
    better_all_extra[[i]] <- rbind(better_all_extra[[i]], data.frame(better_all_extra[[i]][1,1:3] - better_all_extra[[i]][2,1:3],
                                                                model = "diff", author = better_all_extra[[i]]$author[1]))
    better_all_extra[[i]]$is_extra <- "yes"
    temp <- keep_all_auc[[i]]
    temp$is_extra <- "no"
    better_all_extra[[i]] <- rbind(better_all_extra[[i]], temp)
  } else {
    better_all_extra[[i]] <- keep_all_auc[[i]]
    better_all_extra[[i]]$is_extra <- "no"
  }
}

better_all_extra <- do.call("rbind", better_all_extra)
better_all_extra$disease <- convert_names(better_all_extra$author)

better_plot_df <- better_all_extra[better_all_extra$model %in% c("score", "base"),]
ranked_disease <- better_plot_df$disease[better_plot_df$model == "score" & better_plot_df$is_extra == "no"][
  order(better_plot_df$est[better_plot_df$model == "score" & better_plot_df$is_extra == "no"])]
better_plot_df$disease <- factor(better_plot_df$disease, levels = ranked_disease)
better_plot_df$model <- str_to_title(better_plot_df$model)
better_plot_df$is_extra <- str_to_title(better_plot_df$is_extra)
the_plot <- ggplot(better_plot_df, aes(est, disease, color = model, shape = is_extra)) + geom_point(stroke = 1.3) +
  scale_shape_manual(values=c(1, 4)) +
  labs(color = "Model", y = "", x = "AUC", shape = "Includes\nExtra\nCovars.") +
  geom_hline(yintercept = 1:length(unique(better_plot_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("test_plots/meta/other.covar.val.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 8, height = 7)

write_extra <- data.frame("disease" = better_plot_df$disease[better_plot_df$model == "Score" & better_plot_df$is_extra == "No"],
                          "score_yes" = rep(NA, 23),
                          "score_no" = signif(better_plot_df$est[better_plot_df$model == "Score" & better_plot_df$is_extra == "No"], 3),
                          "base_yes" = rep(NA, 23),
                          "base_no" = signif(better_plot_df$est[better_plot_df$model == "Base" & better_plot_df$is_extra == "No"], 3))
write_extra$score_yes[write_extra$disease %in% unique(better_plot_df$disease[better_plot_df$is_extra == "Yes"])] <-
  signif(better_plot_df$est[better_plot_df$model == "Score" & better_plot_df$is_extra == "Yes"], 3)
write_extra$base_yes[write_extra$disease %in% unique(better_plot_df$disease[better_plot_df$is_extra == "Yes"])] <-
  signif(better_plot_df$est[better_plot_df$model == "Base" & better_plot_df$is_extra == "Yes"], 3)
write.table(write_extra, "supp_tables/extra_test.txt", row.names = F, col.names = T, quote = F, sep = '\t')


better_plot_df <- better_all_extra[better_all_extra$model == "diff",]
ranked_disease <- better_plot_df$disease[better_plot_df$is_extra == "no"][order(better_plot_df$est[better_plot_df$is_extra == "no"])]
better_plot_df$disease <- factor(better_plot_df$disease, levels = ranked_disease)
better_plot_df$model <- str_to_title(better_plot_df$model)
better_plot_df$is_extra <- str_to_title(better_plot_df$is_extra)
the_plot <- ggplot(better_plot_df, aes(est, disease, color = is_extra)) + geom_point() +
  labs(y = "", x = "AUC", color = "Includes\nExtra\nCovars.") +
  geom_hline(yintercept = 1:length(unique(better_plot_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("test_plots/meta/other.covar.diff.png"),
       the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 6, height = 5.5)

#stats
better_plot_df <- better_all_extra[better_all_extra$model %in% c("score", "base"),]
temp_df <- better_plot_df[better_plot_df$author %in% names(table(better_plot_df$author)[table(better_plot_df$author) == 4]),]
vals1 <- rep(0, length(unique(temp_df$author)))
vals2 <- rep(0, length(unique(temp_df$author)))
vals3 <- rep(0, length(unique(temp_df$author)))
for(i in 1:length(unique(temp_df$author))){
  auth <- unique(temp_df$author)[i]
  temp_val <- temp_df$est[temp_df$author == auth & temp_df$model == "score"]
  vals1[i] <- temp_val[1] - temp_val[2]
  temp_val <- temp_df$est[temp_df$author == auth & temp_df$is_extra == "yes"]
  vals2[i] <- temp_val[1] - temp_val[2]
  temp_val <- temp_df$est[temp_df$author == auth & temp_df$is_extra == "no"]
  vals3[i] <- temp_val[1] - temp_val[2]
}

mean(vals1)
sd(vals1)/sqrt(length(vals1))
mean(vals3 - vals2)
sd(vals3 - vals2)/sqrt(length(vals1))



################################################################################################
#                                ROC TEST                        ####
#################################################################################################

roc_test_list_basic <- list()
roc_test_list_extra <- list()
for(i in 1:length(best_author)){
  temp <- readRDS(paste0("test_results/", best_author[i], ".roc_data.RDS"))
  roc_test_list_basic[[i]] <- temp[[1]]
  if(!is.null(temp[[3]])){
    roc_test_list_extra[[i]] <- temp[[3]]
  } else {
    roc_test_list_extra[[i]] <- NA
  }
}

basic_roc_pval = data.frame("disease" = convert_names(best_author, author_dict),
                            "roc_pval" = unlist(roc_test_list_basic))
extra_roc_pval = data.frame("disease" = convert_names(best_author, author_dict),
                            "roc_pval" = unlist(roc_test_list_extra))

basic_roc_pval[,2] <- signif(basic_roc_pval[,2], 3)
write.table(basic_roc_pval, "supp_tables/test_roc_pva.txt", sep = "\t", row.names = F, col.names = T, quote = F)

sum(basic_roc_pval$roc_pval < 0.0001)
