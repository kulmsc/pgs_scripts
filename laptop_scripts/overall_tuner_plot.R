library(stringr)
library(ggplot2)
library(cowplot)
library(plyr)
library(viridis)
theme_set(theme_cowplot())



all_best_names <- grep("RDS", list.files("tune_results/", "best"), value = T)
best_author <-  str_split(all_best_names, fixed("."), simplify = T)[,1]
best_author <- best_author[best_author != "Xie"]

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)

convert_names <- function(x, the_dict = method_dict){
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



################################################################################################
#                                 SET UP DATA                                         ####
#################################################################################################

for(i in 1:length(best_author)){

  curr_best <- readRDS(paste0("tune_results/", all_best_names[i]))
  
  all_auc[[i]] <- curr_best[[2]]
  all_auc[[i]] <- all_auc[[i]][!is.na(all_auc[[i]]$auc),]
  all_auc[[i]]$author <- best_author[i]
  all_auc[[i]]$rank <- rank(-1*all_auc[[i]]$auc)
  
  all_conc[[i]] <- curr_best[[1]]
  all_conc[[i]] <- all_conc[[i]][!is.na(all_conc[[i]]$conc),]
  all_conc[[i]]$author <- best_author[i]
  all_conc[[i]]$rank <- rank(-1*all_conc[[i]]$conc)
  
  all_or[[i]] <- curr_best[[4]]
  all_or[[i]] <- all_or[[i]][!is.na(all_or[[i]]$or),]
  all_or[[i]]$author <- best_author[i]
  all_or[[i]]$rank <- rank(-1*all_or[[i]]$or)
  
  all_km[[i]] <- curr_best[[3]]
  all_km[[i]] <- all_km[[i]][!is.na(all_km[[i]]$mean_vals),]
  all_km[[i]]$author <- best_author[i]
  all_km[[i]]$rank <- rank(-1*all_km[[i]]$mean_vals)

}

all_auc <- do.call("rbind", all_auc)
all_conc <- do.call("rbind", all_conc)
all_or <- do.call("rbind", all_or)
all_km <- do.call("rbind", all_km)

# all_auc <- all_auc[!is.na(all_auc$auc),]
# all_or <- all_or[!is.na(all_or$or),]
# all_conc <- all_conc[!is.na(all_conc$conc),]
# all_km <- all_km[!is.na(all_km$mean_vals),]

# all_auc$rank <- abs(all_auc$rank - 15) + 1
# all_conc$rank <- abs(all_conc$rank - 15) + 1
# all_or$rank <- abs(all_km$rank - 15) + 1

# all_auc$method <- convert_names(all_auc$method)
# all_conc$method <- convert_names(all_conc$method)
# all_or$method_name <- convert_names(all_or$method_name)
# all_km$method_name <- convert_names(all_km$method_name)

################################################################################################
#                                 make_function                                    ####
#################################################################################################

make_plots <- function(df, stat_name, plot_name, methods_ordered = order_auc_names){
  
  colnames(df)[colnames(df) == stat_name] <- "stat"

  mean_method_stat <- plyr::daply(df, .(method), function(x) mean(x$stat))
  df$method <- factor(df$method, levels = order_auc_names)
  
  plot_one <- ggplot(df, aes(stat, method)) + geom_boxplot() +
    labs(x = plot_name, y = "") 
    
  
  mean_method_stat <- plyr::daply(df, .(method), function(x) mean(x$rank))
  df$method <- factor(df$method, levels = names(mean_method_stat)[order(mean_method_stat)])
  #df$rank <- abs(df$rank - max(df$rank) - 1)
  
  plot_two <- ggplot(df, aes(rank, method)) + geom_boxplot() +
    labs(x = paste0(plot_name, " Rank"), y = "") 
    

  
  plot_three <- ggplot(df, aes(rank, method, color = stat)) + 
    geom_point(position = position_jitter(width = 0, height = 0.3)) +
    labs(x = paste0(plot_name, " Rank"), y = "", color = plot_name) +
    scale_color_viridis() +
    geom_hline(yintercept = 1:length(unique(df$method)) + 0.5, color = "grey80")
  
  
  df$for_group <- paste0(df$method, "_", df$rank)
  sub_stat <- df[rep(1, length(unique(df$for_group))),]
  sub_stat$domt <- 0
  i <- 1
  for(ugroup in unique(df$for_group)){
    sub_stat$method[i] <- df$method[df$for_group == ugroup][1]
    sub_stat$rank[i] <- df$rank[df$for_group == ugroup][1]
    sub_stat$stat[i] <- mean(df$stat[df$for_group == ugroup], na.rm = T)
    sub_stat$domt[i] <- sum(df$for_group == ugroup)
    i <- i + 1
  }
  
  plot_four <- ggplot(sub_stat, aes(rank, method, color = domt)) + 
    geom_point(size = 2) +
    labs(x = paste0(plot_name, " Rank"), y = "", color = "Number\nDiseases") +
    scale_color_viridis() +
    geom_hline(yintercept = 1:length(unique(df$method)) + 0.5, color = "grey80")

  return(list(plot_one, plot_two, plot_three, plot_four))
}


save_func <- function(plot_list, inds, overall_class, plot_types){
  j <- 1
  for(i in inds){
    ggsave(paste0("tune_plots/meta/", overall_class, ".",  plot_types[j], ".png"), plot_list[[i]], "png", width = 5, height = 5)
    j <- j + 1
  }
}

std <- function(x) sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))

get_table <- function(df, stat_name){
  outdf <- data.frame(matrix(0, nrow = length(unique(df$method)), ncol = 5))
  umeth <- unique(df$method)
  outdf[,1] <- umeth
  df[is.infinite(df[,which(colnames(df) == stat_name)]),which(colnames(df) == stat_name)] <- NA
  for(i in 1:length(umeth)){
    outdf[i,2] <- signif(mean(df[df$method == umeth[i], which(colnames(df) == stat_name)], na.rm = T), 3)
    outdf[i,3] <- signif(std(df[df$method == umeth[i], which(colnames(df) == stat_name)]), 3)
    outdf[i,4] <- signif(mean(df[df$method == umeth[i], which(colnames(df) == "rank")], na.rm = T), 3)
    outdf[i,5] <- signif(std(df[df$method == umeth[i], which(colnames(df) == "rank")]), 3)
  }
  colnames(outdf) <- c("Method", "Stat_Mean", "Stat_SE", "Rank_Mean", "Rank_SE")
  return(outdf)
}

auc_table <- get_table(all_auc, "auc")
write.table(auc_table, "supp_tables/auc_tune.txt", row.names = F, col.names = T, quote = F, sep = "\t")
or_table <- get_table(all_or, "or")
write.table(or_table, "supp_tables/or_tune.txt", row.names = F, col.names = T, quote = F, sep = "\t") 
conc_table <- get_table(all_conc, "conc")
write.table(conc_table, "supp_tables/conc_tune.txt", row.names = F, col.names = T, quote = F, sep = "\t") 
km_table <- get_table(all_km, "mean_vals")
write.table(km_table, "supp_tables/km_tune.txt", row.names = F, col.names = T, quote = F, sep = "\t") 

#written stats
table(all_auc$rank[all_auc$method == "prsCS"])

overall_rank_mean <- cbind((auc_table$Rank_Mean + or_table$Rank_Mean + conc_table$Rank_Mean + km_table$Rank_Mean)/4,
                           as.character(auc_table$Method))

chris_comp <- all_auc[all_auc$author == "Christophersen" & all_auc$method %in% c("prsCS", "WC-2D"),]

################# Quick Correlations #######################3

all_cors <- c(cor(all_auc$rank, all_conc$rank, method = "spearman"),
              cor(all_auc$rank, all_or$rank, method = "spearman"),
              cor(all_auc$rank, all_km$rank, method = "spearman"),
              cor(all_conc$rank, all_or$rank, method = "spearman"),
              cor(all_conc$rank, all_km$rank, method = "spearman"),
              cor(all_or$rank, all_km$rank, method = "spearman"))
all_cors_pval <- c(cor.test(all_auc$rank, all_conc$rank, method = "spearman")$p.value,
              cor.test(all_auc$rank, all_or$rank, method = "spearman")$p.value,
              cor.test(all_auc$rank, all_km$rank, method = "spearman")$p.value,
              cor.test(all_conc$rank, all_or$rank, method = "spearman")$p.value,
              cor.test(all_conc$rank, all_km$rank, method = "spearman")$p.value,
              cor.test(all_or$rank, all_km$rank, method = "spearman")$p.value)


################################################################################################
#                                 AUC                                ####
#################################################################################################
# 0.0057709
for(x in unique(all_auc$author)){
  print(all_auc$auc[all_auc$author == x & all_auc$method == "prsCS"] - all_auc$auc[all_auc$author == x & all_auc$rank == 1])
}

mean_method_stat <- plyr::daply(all_auc, .(method), function(x) mean(x$auc))
order_auc_names <- names(mean_method_stat)[order(mean_method_stat)]

the_plots <- make_plots(all_auc, "auc", "AUC")
save_func(the_plots, c(1, 2, 3, 4), c("abs_auc"), c("val", "rank", "all_shown", "xchosen"))
ggsave("tune_plots/meta/fig1a.auc.xchosen.png", the_plots[[4]] + theme(legend.position = "top") + theme(legend.key.width = unit(1, "cm")),
        "png", width = 4.5, height = 6)

the_plots2 <- make_plots(all_auc, "diff_auc", "AUC Improvement")
save_func(the_plots2, c(1, 2, 3, 4), c("diff_auc"), c("val", "rank", "all_shown", "xchosen"))


################################################################################################
#                                 CONC                                       ####
#################################################################################################


the_plots <- make_plots(all_conc, "conc", "Concordance")
save_func(the_plots, c(1, 2, 4), c("abs_conc"), c("val", "rank", "xchosen"))


the_plots2 <- make_plots(all_conc, "diff_conc", "Concordance Improvement")
save_func(the_plots, c(1, 4), c("diff_conc"), c("val", "xchosen"))



################################################################################################
#                                 OR                                        ####
#################################################################################################

colnames(all_or)[4:5] <- c("score_name", "method")
all_or <- all_or[!is.na(all_or$or) & !is.infinite(all_or$or),]
the_plots <- make_plots(all_or, "or", "Odds Ratio")
save_func(the_plots, c(1, 2, 4), c("abs_or"), c("val", "rank", "xchosen"))


the_plots2 <- make_plots(all_or, "diff_or", "Odds Ratio")
save_func(the_plots, c(1, 4), c("diff_or"), c("val", "xchosen"))

################################################################################################
#                                 KM                                        ####
#################################################################################################

colnames(all_km)[4:5] <- c("score_name", "method")
the_plots <- make_plots(all_km, "mean_vals", "Disease Incidence\nDifference")
save_func(the_plots, c(1, 2, 4), c("abs_km"), c("val", "rank", "xchosen"))

the_plots <- make_plots(all_km, "diff_mean_vals", "Disease Incidence\nDifference")
save_func(the_plots, c(1, 4), c("diff_km"), c("val","xchosen"))

################################################################################################
#                                 ALL TOGETHER                                        ####
#################################################################################################


#do point, boxplots of ranks
all_auc$name_stat <- "AUC"
all_conc$name_stat <- "Concordance"
all_or$name_stat <- "Odds Ratio"
all_km$name_stat <- "Disease Incidence\nDifference"

new_names <- c("score_name", "method", "stat", "rank", "name_stat")
all_toget <- rbind(setNames(all_auc[c("score_name", "method", "auc", "rank", "name_stat")], new_names),
                   setNames(all_conc[c("score_name", "method", "conc", "rank", "name_stat")], new_names),
                   setNames(all_or[c("score_name", "method", "or", "rank", "name_stat")], new_names),
                   setNames(all_km[c("score_name", "method", "mean_vals", "rank", "name_stat")], new_names))
all_toget$score_name <- as.character(all_toget$score_name)
all_toget$method <- as.character(all_toget$method)

summary_toget <- all_toget[rep(1, length(unique(all_toget$method)) * length(unique(all_toget$name_stat))),]
summary_toget$method <- rep(unique(all_toget$method), length(unique(all_toget$name_stat)))
summary_toget$name_stat <- rep(unique(all_toget$name_stat), each = length(unique(all_toget$method)))
summary_toget$spread <- 0

for(i in 1:nrow(summary_toget)){
  umethod <- summary_toget$method[i]
  ustat <- summary_toget$name_stat[i]
  summary_toget$rank[i] <- median(all_toget$rank[all_toget$method == umethod & all_toget$name_stat == ustat])
  summary_toget$spread[i] <- summary_toget$rank[i] - quantile(all_toget$rank[all_toget$method == umethod & all_toget$name_stat == ustat], 0.25)
}

method_ordered <- summary_toget$method[summary_toget$name_stat == "AUC"][
  order(summary_toget$rank[summary_toget$name_stat == "AUC"])]
all_toget$method <- factor(all_toget$method, levels = method_ordered)
summary_toget$method <- factor(summary_toget$method, levels = method_ordered)

#do point, boxplots of scaled values
the_plot <- ggplot(all_toget, aes(rank, method, fill = name_stat)) + geom_boxplot() +
  labs(x = "Rank", fill = "Statistic", y = "") +
  geom_hline(yintercept = 1:length(unique(all_toget$method)) + 0.5)
ggsave("tune_plots/meta/all_together.box.png", the_plot, "png", width = 5, height = 5)


the_plot <- ggplot(summary_toget, aes(rank, method, color = name_stat)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.6)) +
  geom_errorbarh(aes(xmax = rank + spread, xmin = rank - spread, height = 0), 
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.6)) +
  geom_hline(yintercept = 1:length(unique(summary_toget$method)) + 0.5) +
  labs(x = "Rank", y = "", color = "Statistic")
ggsave("tune_plots/meta/all_together.dot.png", the_plot, "png", width = 5, height = 5)


summary_toget$name_stat[summary_toget$name_stat == "Disease Incidence\nDifference"] <- "Disease\nIncidence\nDifference"
the_plot <- ggplot(summary_toget[summary_toget$name_stat != "AUC",], aes(rank, method, color = name_stat)) +
  geom_point(position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.6)) +
  geom_errorbarh(aes(xmax = rank + spread, xmin = rank - spread, height = 0), 
                 position = position_jitterdodge(jitter.width = 0, jitter.height = 0, dodge.width = 0.6)) +
  geom_hline(yintercept = 1:length(unique(summary_toget$method)) + 0.5, color = "grey80") +
  labs(x = "Rank", y = "", color = "Statistic")
ggsave("tune_plots/meta/fig1b.summ.png", the_plot +
         theme(legend.position = "top", plot.margin=unit(c(1,2,1,1), "cm")) +
         guides(color=guide_legend(ncol=3)) +
         theme(legend.text=element_text(size=8)), 
       "png", width = 5, height = 6)
