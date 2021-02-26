library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(pROC)
theme_set(theme_cowplot())

#all plots should look across all traits
#some plots should be looking at the best methods for each trait
#some plots should look at only the single best method for a trait
#should try and replicate the plots made for each individual disease, but now we take the mean in a different way

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
author_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = '\t')
all_authors <- str_split(list.files("nonpred_results/per_score/"), fixed("."), simplify = T)[,1]
all_authors <- unique(all_authors[all_authors != "xie"])
author_dict[,1] <- tolower(author_dict[,1])

convert_names <- function(x, the_dict = author_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}




get_plot_df <- function(rds_name){
    temp_list <- list()
    for(i in 1:length(all_authors)){
      if(file.exists(paste0("nonpred_results/derived_from_per_score/", all_authors[i], ".", rds_name, ".RDS"))){
        temp_list[[i]] <- readRDS(paste0("nonpred_results/derived_from_per_score/", all_authors[i], ".", rds_name, ".RDS"))
        temp_list[[i]]$author <- all_authors[i]
        
        best_method <- read.table(paste0("tune_results/", all_authors[i], ".methods.ss"))
        method_number <- str_split(best_method[,1], fixed("."), simplify = T)[,2]
        method_name <- convert_names(str_split(best_method[,1], fixed("."), simplify = T)[,3], method_dict)
        if(rds_name == "ethnic.plot"){
          temp_list[[i]]$neat_score_names <- paste0(temp_list[[i]]$method, "-", str_split(temp_list[[i]]$name, fixed("."), simplify = T)[,2])
        }
        if(rds_name == "sibs_res.plot"){
          temp_list[[i]]$neat_score_names <- paste0(temp_list[[i]]$method, "-",
                                                    str_split(row.names(temp_list[[i]]), fixed("."), simplify = T)[,2])
        }
        if(rds_name == "ethnic.tstat" | rds_name == "ethnic.pval"){
          temp_list[[i]] <- melt(temp_list[[i]], id.vars = c("eth1", "eth2"))
          temp_list[[i]] <- temp_list[[i]][temp_list[[i]]$variable != "author",]
          temp_list[[i]]$method_number <- unlist(lapply(strsplit(as.character(temp_list[[i]]$variable), ".", fixed = T), function(x) x[2]))
          temp_list[[i]]$method_name <- unlist(lapply(strsplit(as.character(temp_list[[i]]$variable), ".", fixed = T), function(x) x[3]))
          temp_list[[i]]$method_name <- convert_names(temp_list[[i]]$method_name, method_dict)
          temp_list[[i]]$neat_score_names <- paste0(temp_list[[i]]$method_name, "-", temp_list[[i]]$method_number)
          temp_list[[i]]$author <- all_authors[i]
        }
        if(rds_name != "ethnic.plot" & rds_name != "ethnic.tstat" & rds_name != "ethnic.pval"){
          if("neat_method_names" %in% colnames(temp_list[[i]]) | "neat_score_names" %in% colnames(temp_list[[i]])){
            temp_list[[i]] <- temp_list[[i]][temp_list[[i]]$neat_score_names %in% paste0(method_name, "-", method_number),]
          } else {
            temp_list[[i]] <- temp_list[[i]][temp_list[[i]]$score %in% paste0(method_name, "-", method_number),]
          }
        }
        
        best_method <- read.table(paste0("tune_results/", all_authors[i], ".best.ss"))
        method_number <- str_split(best_method[,2], fixed("."), simplify = T)[,2]
        method_name <- convert_names(str_split(best_method[,2], fixed("."), simplify = T)[,3], method_dict)
        temp_list[[i]]$best <- 0 
        if("neat_method_names" %in% colnames(temp_list[[i]]) | "neat_score_names" %in% colnames(temp_list[[i]])){
          temp_list[[i]]$best[temp_list[[i]]$neat_score_names %in% paste0(method_name[3], "-", method_number[3])] <- 1
        } else {
          temp_list[[i]]$best[temp_list[[i]]$score %in% paste0(method_name[3], "-", method_number[3])] <- 1
        }
        if(sum(temp_list[[i]]$best) == 0){
          temp_list[[i]]$best[temp_list[[i]]$method == method_name[3]] <- 1
        }
        if(sum(temp_list[[i]]$best) == 0){
          temp_list[[i]]$best[temp_list[[i]]$method == "prsCS"] <- 1
        }
        
      }
    }
    temp_list <- do.call("rbind", temp_list)
    temp_list$disease <- convert_names(temp_list$author, author_dict)
    
    return(temp_list)
}


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# get base AUCS
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
all_files <- list.files("test_results/",  "res")
all_files <- all_files[-1]
all_base_auc <- rep(0, 23)
for(i in 1:length(all_files)){
  print(i)
  res <- readRDS(paste0("test_results/", all_files[i]))
  all_base_auc[i] <- res[["base"]][["auc"]][[2]][2]
}
all_base_auc <- data.frame("disease" = convert_names(tolower(str_split(all_files, fixed("_"), simplify = T)[,1]), author_dict),
                           "auc" = all_base_auc)
all_base_auc$disease <- as.character(all_base_auc$disease)

########################################################################################################################
#                                            MOD   SET    STUFF
#########################################################################################################################

###########################################################
#              SCORE SIZES                              ##
###########################################################

score_size_df <- get_plot_df("score.len")
score_size_df$log_len <- log10(score_size_df$all_len)

mean_size <- plyr::daply(score_size_df, .(disease), function(x) median(x$log_len))
score_size_df$disease <- factor(score_size_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_size_df, aes(all_len, disease)) + geom_boxplot() +
  labs(x = "log10(Score Length)", y = "") +
  scale_x_continuous(trans='log10') +
  geom_hline(yintercept = 1:length(unique(score_size_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_size.1.png", the_plot, height = 5, width = 6.5)


temp_df <- score_size_df[score_size_df$best == 1,]
temp_df <- temp_df[!duplicated(temp_df$author),]
temp_df$disease <- factor(temp_df$disease, levels = temp_df$disease[order(temp_df$all_len)])
the_plot <- ggplot(temp_df, aes(all_len, disease)) + geom_point() +
  scale_x_continuous(trans='log10') +
  labs(x = "log10(Score Length)", y = "")  +
  geom_hline(yintercept = 1:length(unique(score_size_df$disease)) + 0.5, color = "grey80") +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
plot(the_plot)
ggsave("meta_plots/meta/score_size.2.png", the_plot, height = 5, width = 6.5)
write.table(temp_df[,c(6,3,7)], "supp_tables/score_len.txt", row.names = F, col.names = T, quote = F, sep = "\t")


mean_size <- plyr::daply(score_size_df, .(neat_method_names), function(x) median(x$log_len))
score_size_df$disease <- factor(score_size_df$neat_method_names, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_size_df, aes(all_len, neat_method_names)) + geom_boxplot() +
  labs(x = "log10(Score Length)", y = "") +
  scale_x_continuous(trans='log10') +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))
plot(the_plot)
ggsave("meta_plots/meta/score_size.3.png", the_plot, height = 5, width = 6.5)




###########################################################
#              SCORE MEAN EFFECT                              ##
###########################################################

score_quant_df <- get_plot_df("score.mean_quant")

mean_size <- plyr::daply(score_quant_df[score_quant_df$variable == "Q4",], .(disease), function(x) median(x$value, na.rm = T))
score_quant_df$disease <- factor(score_quant_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_quant_df, aes(value, disease, color = variable)) + geom_boxplot() +
  labs(x = "Mean Effect", y = "", color = "Group") +
  geom_hline(yintercept = 1:length(unique(score_quant_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_mean.1.png", the_plot, height = 5, width = 6.5)

the_plot <- ggplot(score_quant_df[score_quant_df$best == 1,], aes(value, disease, color = variable)) + geom_point() +
  labs(x = "Mean Effect", y = "", color = "Group") +
  geom_hline(yintercept = 1:length(unique(score_quant_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_mean.2.png", the_plot, height = 5, width = 6.5)
temp_df <- score_quant_df[score_quant_df$best == 1,]
temp_df <- data.frame("Disease" = temp_df$disease[temp_df$variable == "Q1"],
                      "Q1" = signif(temp_df$value[temp_df$variable == "Q1"], 3),
                      "Q2" = signif(temp_df$value[temp_df$variable == "Q2"], 3),
                      "Q3" = signif(temp_df$value[temp_df$variable == "Q3"], 3),
                      "Q4" = signif(temp_df$value[temp_df$variable == "Q4"], 3))
write.table(temp_df, "supp_tables/score_lengroup.txt", row.names = F, col.names = T, quote = F, sep = "\t")

mean_size <- plyr::daply(score_quant_df[score_quant_df$variable == "Q4",], .(method), function(x) median(x$value, na.rm = T))
score_quant_df$method <- factor(score_quant_df$method, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_quant_df, aes(value, method, color = variable)) + geom_boxplot() +
  labs(x = "Mean Effect", y = "", color = "Group") +
  geom_hline(yintercept = 1:length(unique(score_quant_df$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_mean.3.png", the_plot, height = 5, width = 6.5)



###########################################################
#              SCORE MEAN LENGTH                              ##
###########################################################

score_len_df <- get_plot_df("score.len_quant")


mean_size <- plyr::daply(score_len_df[score_len_df$variable == "Q4",], .(disease), function(x) median(x$value, na.rm = T))
score_len_df$disease <- factor(score_len_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_len_df, aes(value, disease, color = variable)) + geom_boxplot() +
  labs(x = "Mean Size", y = "", color = "Group") +
  geom_hline(yintercept = 1:length(unique(score_len_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_len.1.png", the_plot, height = 5, width = 6.5)

plot_df <- score_len_df[score_len_df$best == 1,]
plot_df$disease <- as.character(plot_df$disease)
plot_df$disease <- factor(plot_df$disease, levels = plot_df$disease[plot_df$variable == "Q1"][order(plot_df$value[plot_df$variable == "Q1"])])
the_plot <- ggplot(plot_df, aes(value, disease, color = variable)) + geom_point() +
  labs(x = "Mean Size", y = "", color = "Equal\nEffect\nSize\nGroup") +
  geom_hline(yintercept = 1:length(unique(score_len_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_len.2.png", the_plot, height = 5, width = 6.5)

mean_size <- plyr::daply(score_len_df[score_len_df$variable == "Q4",], .(method), function(x) median(x$value, na.rm = T))
score_len_df$method <- factor(score_len_df$method, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_len_df, aes(value, method, color = variable)) + geom_boxplot() +
  labs(x = "Mean Size", y = "", color = "Group") +
  geom_hline(yintercept = 1:length(unique(score_len_df$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/score_len.3.png", the_plot, height = 5, width = 6.5)




########################################################################################################################
#                                        PHENO    DEF
#########################################################################################################################

score_pheno_def_df <- get_plot_df("pheno_def.auc")
#temp <- paste0(score_pheno_def_df$value, score_pheno_def_df$method, score_pheno_def_df$author, score_pheno_def_df$variable)
#score_pheno_def_df <- score_pheno_def_df[!duplicated(temp),]


score_pheno_def_df$diff <- 0 
score_pheno_def_df$variable <- as.character(score_pheno_def_df$variable)
score_pheno_def_df$variable[score_pheno_def_df$variable == "ICD.or.Self.Rep."] <- "ICD or\nSelf Rep."
for(disease in unique(score_pheno_def_df$disease)){
  for(method in unique(score_pheno_def_df$method)){
    for(vari in c("ICD", "Self\nReported", "Any", "Double\nReported")){
      score_pheno_def_df$diff[score_pheno_def_df$disease == disease & 
                                score_pheno_def_df$method == method & score_pheno_def_df$variable == vari] <-
        score_pheno_def_df$value[score_pheno_def_df$disease == disease & 
                                  score_pheno_def_df$method == method & score_pheno_def_df$variable == vari] - 
        score_pheno_def_df$value[score_pheno_def_df$disease == disease & 
                                  score_pheno_def_df$method == method & score_pheno_def_df$variable == "ICD or\nSelf Rep."]
    }
  }
}

score_pheno_def_df <- score_pheno_def_df[score_pheno_def_df$variable != "ICD or\nSelf Rep.",]

mean_size <- plyr::daply(score_pheno_def_df, .(method), function(x) max(x$value, na.rm = T) - min(x$value, na.rm = T))
score_pheno_def_df$method <- factor(score_pheno_def_df$method, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_pheno_def_df, aes(diff, method, color = variable)) + geom_boxplot() +
  labs(x = "AUC Difference", y = "", color = "Phenotyping\nMethod") +
  geom_hline(yintercept = 1:length(unique(score_pheno_def_df$method)) + 0.5, color = "grey80")
plot(the_plot)
 ggsave("meta_plots/meta/pheno_auc.1.png", the_plot, height = 6, width = 7.5)


mean_size <- plyr::daply(score_pheno_def_df, .(disease), function(x) max(x$value, na.rm = T) - min(x$value, na.rm = T))
score_pheno_def_df$disease <- factor(score_pheno_def_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(score_pheno_def_df, aes(diff, disease, color = variable)) + geom_boxplot() +
  labs(x = "AUC Difference", y = "", color = "Phenotyping\nMethod") +
  geom_hline(yintercept = 1:length(unique(score_pheno_def_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/pheno_auc.2.png", the_plot, height = 6, width = 7.5)

the_plot <- ggplot(score_pheno_def_df, aes(diff, variable)) + geom_boxplot() +
  labs(x = "AUC Difference", y = "")
plot(the_plot)
ggsave("meta_plots/meta/pheno_auc.3.png", the_plot, height = 6, width = 7.5)

temp_df <- score_pheno_def_df[score_pheno_def_df$best == 1,]
temp_df <- temp_df[,c(7,3,4,8)]
temp_df$diff <- signif(temp_df$diff, 3)
temp_df <- data.frame("Disease" = temp_df$disease[temp_df$variable == "ICD"],
                      "AUC ICD" =  temp_df$diff[temp_df$variable == "ICD"],
                      "AUC Self-Reported" = temp_df$diff[temp_df$variable == "Self.Reported"],
                      "AUC Any" = temp_df$diff[temp_df$variable == "Any"],
                      "AUC Double-Reported" = temp_df$diff[temp_df$variable == "Double.Reported"])
write.table(temp_df, "supp_tables/pheno_def.txt", row.names = F, col.names = T, quote = F, sep = "\t")



########################################################################################################################
#                                      SEX     DIFF
#########################################################################################################################

sex_auc_df <- get_plot_df("sex.auc")

mean_size <- plyr::daply(sex_auc_df, .(methods_name), function(x) median(x$male_auc, na.rm = T))
sex_auc_df$methods_name <- factor(sex_auc_df$methods_name, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(sex_auc_df, aes(male_auc, methods_name, color = sex)) + geom_boxplot() +
  labs(x = "AUC", y = "", color = "Sex") +
  geom_hline(yintercept = 1:length(unique(sex_auc_df$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/sex_auc.1.png", the_plot, height = 6, width = 7.5)

mean_size <- plyr::daply(sex_auc_df, .(disease), function(x) median(x$male_auc, na.rm = T))
sex_auc_df$disease <- factor(sex_auc_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(sex_auc_df, aes(male_auc, disease, color = sex)) + geom_boxplot() +
  labs(x = "AUC", y = "", color = "Sex") +
  geom_hline(yintercept = 1:length(unique(sex_auc_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/sex_auc.2.png", the_plot, height = 6, width = 7.5)

the_plot <- ggplot(sex_auc_df[sex_auc_df$best == 1,], aes(male_auc, disease, color = sex)) + geom_point() +
  geom_errorbarh(aes(xmin = male_lo, xmax = male_hi), height = 0) +
  labs(x = "AUC", y = "", color = "Sex") +
  geom_hline(yintercept = 1:length(unique(sex_auc_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/sex_auc.3.png", the_plot, height = 6, width = 6)

for_supp <- sex_auc_df[sex_auc_df$best == 1, c(9,1,2,3,6)]
for_supp[,2] <- signif(for_supp[,2], 3)
for_supp[,3] <- signif(for_supp[,3], 3)
for_supp[,4] <- signif(for_supp[,4], 3)
for_supp <- for_supp[order(for_supp$disease, decreasing = T),]
for_supp <- data.frame("Disease" = for_supp$disease[for_supp$sex == "Male"], 
                       "Male AUC" = for_supp$male_auc[for_supp$sex == "Male"],
                       "Male AUC CI" = abs(for_supp$male_auc[for_supp$sex == "Male"] - for_supp$male_hi[for_supp$sex == "Male"]),
                       "Female AUC" = for_supp$male_auc[for_supp$sex == "Female"],
                       "Female AUC CI" = abs(for_supp$male_auc[for_supp$sex == "Female"] - for_supp$male_hi[for_supp$sex == "Female"]))
for_supp <- for_supp[!is.na(for_supp$Male.AUC),]
write.table(for_supp, "supp_tables/sex_diff.txt", row.names = F, col.names = T, quote = F, sep = "\t")


comp_to_base <- function(strat_df, base_df){
  for(d in unique(strat_df$disease)){
    
  }
}

########################################################################################################################
#                                      ETHNICITY DIFF
#########################################################################################################################



ethnic_auc_df <- get_plot_df("ethnic.plot")

ethnic_auc_df$diff <- 0 
for(disease in unique(ethnic_auc_df$disease)){
  for(method in unique(ethnic_auc_df$method)){
    for(vari in c("African", "Asian", "European")){
        ethnic_auc_df$diff[ethnic_auc_df$disease == disease & 
                                ethnic_auc_df$method == method & ethnic_auc_df$ethnic == vari] <-
        ethnic_auc_df$mean[ethnic_auc_df$disease == disease & 
                                   ethnic_auc_df$method == method & ethnic_auc_df$ethnic == vari] - 
        ethnic_auc_df$mean[ethnic_auc_df$disease == disease & 
                                   ethnic_auc_df$method == method & ethnic_auc_df$ethnic == "British"]
    }
  }
}

ethnic_auc_df <- ethnic_auc_df[ethnic_auc_df$ethnic != "British",]

the_plot <- ggplot(ethnic_auc_df, aes(diff, method, color = ethnic)) + geom_boxplot() +
  labs(x = "Difference from UK Score", y = "", color = "Population\nGroup") +
  geom_hline(yintercept = 1:length(unique(ethnic_auc_df$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/ethnic_auc.1.png", the_plot, height = 6, width = 6)


the_plot <- ggplot(ethnic_auc_df, aes(diff, disease, color = ethnic)) + geom_boxplot() +
  labs(x = "Difference from UK Score", y = "", color = "Population\nGroup") +
  geom_hline(yintercept = 1:length(unique(ethnic_auc_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/ethnic_auc.2.png", the_plot, height = 6, width = 6)

mean_size <- plyr::daply(ethnic_auc_df[ethnic_auc_df$best == 1,], .(disease), function(x) var(x$diff))
ethnic_auc_df$disease <- factor(ethnic_auc_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(ethnic_auc_df[ethnic_auc_df$best == 1,], aes(diff, disease, color = ethnic)) + geom_point() +
  labs(x = "Difference from UK Score", y = "", color = "Population\nGroup") +
  geom_hline(yintercept = 1:length(unique(ethnic_auc_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/ethnic_auc.3.png", the_plot, height = 6, width = 6)

temp_df <- ethnic_auc_df[ethnic_auc_df$best == 1,]
temp_df <- temp_df[,c(9, 1, 2, 10, 3)]
temp_df$diff <- signif(temp_df$diff, 3)
temp_df <- data.frame("Disease" = temp_df$disease[temp_df$ethnic == "European"],
                      "European" = temp_df$diff[temp_df$ethnic == "European"],
                      "African" = temp_df$diff[temp_df$ethnic == "African"],
                      "Asian" = temp_df$diff[temp_df$ethnic == "Asian"])
write.table(temp_df, "supp_tables/ethnic_dists.txt", row.names = F, col.names = T, quote = F, sep = "\t")


temp_df$diff[temp_df$disease == "Ovarian Cancer"]/temp_df$mean[temp_df$disease == "Ovarian Cancer"]
temp_df$diff[temp_df$disease == "NAFLD"]/temp_df$mean[temp_df$disease == "NAFLD"]
mean(abs(temp_df$diff[temp_df$ethnic == "European"]/temp_df$mean[temp_df$ethnic == "European"]))

#############
ethnic_tstat_df <- get_plot_df("ethnic.tstat")
ethnic_tstat_df$comp <- paste0(ethnic_tstat_df$eth1, "-", ethnic_tstat_df$eth2)
ethnic_tstat_df$value <- as.numeric(ethnic_tstat_df$value)

the_plot <- ggplot(ethnic_tstat_df[ethnic_tstat_df$best == 1, ], aes(value, disease, color = comp)) +
  geom_vline(xintercept = 0) + geom_point() +
  labs(x = "T - Statistic", y = "", color = "Pop. Group\nComparison") +
  geom_hline(yintercept = 1:length(unique(ethnic_tstat_df$disease)) + 0.5, color = "grey80") 
plot(the_plot)
ggsave("meta_plots/meta/ethnic_auc.4.png", the_plot, height = 6, width = 6)

ethnic_pval_df <- get_plot_df("ethnic.pval")
ethnic_pval_df$comp <- paste0(ethnic_pval_df$eth1, "-", ethnic_pval_df$eth2)
ethnic_pval_df$value <- as.numeric(ethnic_pval_df$value)
ethnic_pval_df$value[ethnic_pval_df$value < 1e-300] <- 1e-300

the_plot <- ggplot(ethnic_pval_df[ethnic_pval_df$best == 1, ], aes(-log10(value), disease, color = comp)) +
  geom_vline(xintercept = 0) + geom_point() +
  labs(x = "p - Value", y = "", color = "Pop. Group\nComparison") +
  geom_hline(yintercept = 1:length(unique(ethnic_pval_df$disease)) + 0.5, color = "grey80") 
plot(the_plot)
ggsave("meta_plots/meta/ethnic_auc.5.png", the_plot, height = 6, width = 6)

temp_df <- ethnic_pval_df[ethnic_pval_df$best == 1, ]
temp_df <- temp_df[,c(10, 1, 2, 4)]
temp_df$eth1 <- convert_names(temp_df$eth1, cbind(c("asian", "african", "euro", "brit"), c("Asian", "African", "European", "British")))
temp_df$eth2 <- convert_names(temp_df$eth2, cbind(c("asian", "african", "euro", "brit"), c("Asian", "African", "European", "British")))
temp_df$value <- signif(temp_df$value, 3)
temp_df <- data.frame("Disease" = temp_df$disease[temp_df$eth1 == "British" & temp_df$eth2 == "European"],
                      "British - European" = temp_df$value[temp_df$eth1 == "British" & temp_df$eth2 == "European"],
                      "British - Asian" = temp_df$value[temp_df$eth1 == "British" & temp_df$eth2 == "Asian"],
                      "British - African" = temp_df$value[temp_df$eth1 == "British" & temp_df$eth2 == "African"],
                      "European - African" = temp_df$value[temp_df$eth1 == "European" & temp_df$eth2 == "African"],
                      "European - Asian" = temp_df$value[temp_df$eth1 == "European" & temp_df$eth2 == "Asian"],
                      "African - Asian" = temp_df$value[temp_df$eth1 == "African" & temp_df$eth2 == "Asian"])
write.table(temp_df, "supp_tables/ethnic_pval.txt", row.names = F, col.names = T, quote = F, sep = "\t")



########################################################################################################################
#                                      SIB    CONC
#########################################################################################################################


sibs_df <- get_plot_df("sibs_res.plot")

sibs_df$diff <- 0 
for(disease in unique(sibs_df$disease)){
  for(method in unique(sibs_df$method)){
    sibs_df$diff[sibs_df$disease == disease & sibs_df$method == method] <- 
      sibs_df$conc[sibs_df$disease == disease & sibs_df$method == method & sibs_df$type == "Sibling"] - 
      sibs_df$conc[sibs_df$disease == disease & sibs_df$method == method & sibs_df$type != "Sibling"]
  }
}

sub_sibs_df <- sibs_df[sibs_df$type == "Sibling",]

mean_size <- plyr::daply(sub_sibs_df, .(disease), function(x) median(x$diff))
sub_sibs_df$disease <- factor(sub_sibs_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(sub_sibs_df, aes(diff, disease)) + geom_boxplot() +
  labs(x = "Sib/Non-Sib Conc. Difference", y = "") + geom_vline(xintercept = 0, color = "black")
plot(the_plot)
ggsave("meta_plots/meta/sibs.conc.1.png", the_plot, height = 6, width = 6)

mean_size <- plyr::daply(sub_sibs_df, .(method), function(x) median(x$diff))
sub_sibs_df$method <- factor(sub_sibs_df$method, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(sub_sibs_df, aes(diff, method)) + geom_boxplot() + geom_vline(xintercept = 0, color = "black")
plot(the_plot)
ggsave("meta_plots/meta/sibs.conc.2.png", the_plot, height = 6, width = 6)

sibs_df$disease <- factor(sibs_df$disease, levels = unique(sibs_df$disease)[order(sibs_df$conc[sibs_df$best == 1 & sibs_df$type == "Sibling"])])
sibs_df <- sibs_df[!is.na(sibs_df$disease),]
the_plot <- ggplot(sibs_df[sibs_df$best == 1,], aes(conc, disease, color = type)) + geom_point() +
  geom_hline(yintercept = 1:length(unique(sibs_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/sibs.conc.3.png", the_plot, height = 6, width = 6)

for_supp <- sibs_df[sibs_df$best == 1,]
for_supp[order(abs(for_supp$diff)),]
for_supp <- for_supp[,c(9,1,2,3,5)]
for_supp[,2] <- signif(for_supp[,2], 3)
for_supp[,3] <- signif(for_supp[,3], 3)
for_supp[,4] <- signif(for_supp[,4], 3)
for_supp <- data.frame("Disease" = for_supp$disease[for_supp$type == "Sibling"],
                       "Sib. Conc." = for_supp$conc[for_supp$type == "Sibling"],
                       "Sib. CI Conc." = for_supp$conc[for_supp$type == "Sibling"] - for_supp$lo[for_supp$type == "Sibling"],
                       "Non-Sib. Conc." = for_supp$conc[for_supp$type == "Non-Sib"],
                       "Non-Sib. CI Conc." = for_supp$conc[for_supp$type == "Non-Sib"] - for_supp$lo[for_supp$type == "Non-Sib"])
write.table(for_supp, "supp_tables/sib_diff.txt", row.names = F, col.names = T, quote = F, sep = "\t")


########################################################################################################################
#                                      AGE     DIFF
#########################################################################################################################

age_auc_df <- get_plot_df("agediff.auc")

mean_size <- plyr::daply(age_auc_df, .(methods_name), function(x) median(x$old_auc, na.rm = T))
age_auc_df$methods_name <- factor(age_auc_df$methods_name, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(age_auc_df, aes(old_auc, methods_name, color = sex)) + geom_boxplot() +
  labs(x = "AUC", y = "", color = "Age\nGroup") +
  geom_hline(yintercept = 1:length(unique(age_auc_df$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/age_auc.1.png", the_plot, height = 6, width = 7.5)

mean_size <- plyr::daply(age_auc_df, .(disease), function(x) median(x$old_auc, na.rm = T))
age_auc_df$disease <- factor(age_auc_df$disease, levels = names(mean_size)[order(mean_size)])
the_plot <- ggplot(age_auc_df, aes(old_auc, disease, color = sex)) + geom_boxplot() +
  labs(x = "AUC", y = "", color = "Age\nGroup") +
  geom_hline(yintercept = 1:length(unique(age_auc_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/age_auc.2.png", the_plot, height = 6, width = 7.5)

the_plot <- ggplot(age_auc_df[age_auc_df$best == 1,], aes(old_auc, disease, color = sex)) + geom_point() +
  geom_errorbarh(aes(xmin = old_lo, xmax = old_hi), height = 0) +
  labs(x = "AUC", y = "", color = "Sex") +
  geom_hline(yintercept = 1:length(unique(age_auc_df$disease)) + 0.5, color = "grey80")
plot(the_plot)
ggsave("meta_plots/meta/age_auc.3.png", the_plot, height = 6, width = 6)

for_supp <- age_auc_df[age_auc_df$best == 1,c(9,1,2,3,6)]
for_supp[,2] <- signif(for_supp[,2], 3)
for_supp[,3] <- signif(for_supp[,3], 3)
for_supp[,4] <- signif(for_supp[,4], 3)
for_supp <- data.frame("Disease" = for_supp$disease[for_supp$sex == "Young"],
                       "AUC Young" = for_supp$old_auc[for_supp$sex == "Young"],
                       "AUC CI Young" = for_supp$old_auc[for_supp$sex == "Young"] - for_supp$old_lo[for_supp$sex == "Young"],
                       "AUC Old" = for_supp$old_auc[for_supp$sex == "Old"],
                       "AUC CI Old" = for_supp$old_auc[for_supp$sex == "Old"] - for_supp$old_lo[for_supp$sex == "Old"])
write.table(for_supp, "supp_tables/age_diff.txt", row.names = F, col.names = T, quote = F, sep = "\t")





########################################################################################################################
#                                      MANY     DIFF
#########################################################################################################################

vars_examined <- c("Time at\nCurrent\nAddress", "Income", "Number in\nHouse", "Year of\nEducation",
                   "Census:\nMedian Age", "Census:\nUnemployment", "Census:\nHealth", "Census:\nPop. Den.")
vars_simple <- c("Time at Address", "Income", "Number in House", "Age at Education",
                 "Census: Median Age", "Census: Unemployment", "Census: Health", "Census: Pop. Den.")
var_splits <- list(c("1-20 years", "20+ years"), c("< 40,000", "> 40,000"), c("1-2 person", "3+ persons"),
                   c("1-19 years", "20+ years"), c("> 42 years", "< 42 years"), c("> 38", "< 38"),
                   c("> 719", "< 719"), c("> 32 persons\n/hectare", "< 32 persons\n/hectare"))

holder_many <- rep(list(list()), 8)

j <- 1
for(author in all_authors){
  author_many <- readRDS(paste0("nonpred_results/derived_from_per_score/", author, ".many.auc.RDS"))
  

  for(i in 1:8){
    holder_many[[i]][[j]] <- author_many[[i]]
    holder_many[[i]][[j]]$disease <- convert_names(author)
    holder_many[[i]][[j]]$best <- 0
    
    best_method <- read.table(paste0("tune_results/", author, ".best.ss"))
    method_number <- str_split(best_method[,2], fixed("."), simplify = T)[,2]
    method_name <- convert_names(str_split(best_method[,2], fixed("."), simplify = T)[,3], method_dict)
    holder_many[[i]][[j]]$best[holder_many[[i]][[j]]$method == method_name[3]] <- 1
  }
  
  j <- j + 1
}

saved_low_p <- list()
saved_interesting <- list()

for(i in 1:8){
  plot_df <- do.call("rbind", holder_many[[i]])
  
  saved_low_p[[i]] <- plot_df[plot_df$pval < 0.05 & plot_df$best == 1,]
  if(sum(plot_df$pval < 0.05 & plot_df$best == 1) > 0){
    saved_low_p[[i]]$var_ind <- i
  }
  
  if(i %in% 1:3){
    saved_interesting[[i]] <- plot_df[plot_df$best == 1,]
    saved_interesting[[i]]$var_ind <- i
  }
  
  mean_size <- plyr::daply(plot_df, .(methods_name), function(x) median(x$upper_auc, na.rm = T))
  plot_df$methods_name <- factor(plot_df$methods_name, levels = names(mean_size)[order(mean_size)])
  the_plot <- ggplot(plot_df, aes(upper_auc, methods_name, color = var)) + geom_boxplot() +
    labs(color = vars_examined[i], x = "AUC", y = "") +
    scale_color_discrete(labels = var_splits[[i]]) +
    geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80")
  plot(the_plot)
  ggsave(paste0("meta_plots/meta/many.", i, ".1.png"), the_plot, height = 6, width = 6)
  
  mean_size <- plyr::daply(plot_df, .(disease), function(x) median(x$upper_auc, na.rm = T))
  plot_df$disease <- factor(plot_df$disease, levels = names(mean_size)[order(mean_size)])
  the_plot <- ggplot(plot_df, aes(upper_auc, disease, color = var)) + geom_boxplot() +
    labs(color = vars_examined[i], x = "AUC", y = "") +
    scale_color_discrete(labels = var_splits[[i]]) +
    geom_hline(yintercept = 1:length(unique(plot_df$disease)) + 0.5, color = "grey80")
  plot(the_plot)
  ggsave(paste0("meta_plots/meta/many.", i, ".2.png"), the_plot, height = 6, width = 6)
  
  plot_df$disease <- factor(plot_df$disease, levels = plot_df$disease[plot_df$best == 1 & plot_df$var == "young"][
    order(plot_df$upper_auc[plot_df$best == 1 & plot_df$var == "young"])])
  the_plot <- ggplot(plot_df[plot_df$best == 1,], aes(upper_auc, disease, color = var)) + geom_point() +
    labs(color = vars_examined[i], x = "AUC", y = "") +
    scale_color_discrete(labels = var_splits[[i]]) +
    geom_errorbarh(aes(xmin = upper_lo, xmax = upper_hi), height = 0) +
    geom_hline(yintercept = 1:length(unique(plot_df$disease)) + 0.5, color = "grey80")
  plot(the_plot)
  ggsave(paste0("meta_plots/meta/many.", i, ".3.png"), the_plot, height = 6, width = 6)
  
  
  
  for_supp <- plot_df[plot_df$best == 1,c(8,1,2,3,6)]
  for_supp[,2] <- signif(for_supp[,2], 3)
  for_supp[,3] <- signif(for_supp[,3], 3)
  for_supp[,4] <- signif(for_supp[,4], 3)
  for_supp <- data.frame("Disease" = for_supp$disease[for_supp$var == "young"],
                         "AUC Young" = for_supp$upper_auc[for_supp$var == "young"],
                         "AUC CI Young" = signif(for_supp$upper_auc[for_supp$var == "young"] - for_supp$upper_lo[for_supp$var == "young"], 3),
                         "AUC Old" = for_supp$upper_auc[for_supp$var == "old"],
                         "AUC CI Old" = signif(for_supp$upper_auc[for_supp$var == "old"] - for_supp$upper_lo[for_supp$var == "old"], 3))
  colnames(for_supp) <- c("Disease", paste0("AUC - ", vars_simple[i], "Low"), paste0("AUC CI - ", vars_simple[i], "Low"),
                          paste0("AUC - ", vars_simple[i], "Hi"), paste0("AUC CI - ", vars_simple[i], "Hi"))
  colnames(for_supp) <- c("Disease", "AUC - Low Group", "AUC CI - Low Group", "AUC - Hi Group", "AUC CI - Hi Group")
  write.table(for_supp, paste0("supp_tables/many.", i, ".txt"), row.names = F, col.names = T, quote = F, sep = "\t")
}


df = do.call("rbind", saved_low_p)
std_diseases <- c("Heart Failure", "Gout", "Prostate Cancer", "CAD", "A. Fib.")
other_diseases <- unique(df$disease[!(df$disease %in% std_diseases)])
nrow(df[df$disease %in% other_diseases & df$var_ind %in% 1:3,])/2
#varind = 1, depression
#varind = 2, depression
#varind = 3, lupus

df = do.call("rbind", saved_interesting)

mean(abs(df[df$disease %in% std_diseases & df$var == "young",2] - df[df$disease %in% std_diseases & df$var != "young",2]))
sd(abs(df[df$disease %in% std_diseases & df$var == "young",2] - df[df$disease %in% std_diseases & df$var != "young",2]))/
  sqrt(sum(df$disease %in% std_diseases & df$var == "young"))

mean(abs(df[df$disease %in% other_diseases & df$var == "young",2] - df[df$disease %in% other_diseases & df$var != "young",2]))
sd(abs(df[df$disease %in% other_diseases & df$var == "young",2] - df[df$disease %in% other_diseases & df$var != "young",2]))/
  sqrt(sum(df$disease %in% other_diseases & df$var == "young"))

########################################################################################################################
#                                      ALL TOGETHER
#########################################################################################################################


use_disease <- readRDS("local_info/diseas.23.RDS")
#MS ALS

score_pheno_def_df <- score_pheno_def_df[score_pheno_def_df$best == 1,]
score_pheno_def_df$disease <- factor(as.character(score_pheno_def_df$disease), levels = rev(sort(use_disease)))
conv_dict <- cbind(c("ICD", "Self.Reported", "x", "Any", "Double.Reported"),
                   c("ICD", "Self\nReported", "ICD or\nSelf Rep.", "Any", "Double\nReported"))
score_pheno_def_df$variable <- convert_names(score_pheno_def_df$variable, conv_dict)

sex_auc_df <- sex_auc_df[sex_auc_df$best == 1,]
add_to <- use_disease[!(use_disease %in% sex_auc_df$disease)]
temp_list <- list()
for(i in 1:length(add_to)){
  temp_list[[i]] <- sex_auc_df[1:2,]
  temp_list[[i]]$disease <- add_to[i]
  temp_list[[i]][1,1:3] <- NA
  temp_list[[i]][2,1:3] <- NA
}

many_aucs <- rep(list(list()), 8)
for(i in 1:length(all_authors)){
  many_raw <- readRDS(paste0("nonpred_results/derived_from_per_score/", all_authors[i], ".many.auc.RDS"))
  best_method <- read.table(paste0("tune_results/", all_authors[i], ".best.ss"))
  best_method <- convert_names(strsplit(as.character(best_method[3,2]), ".", fixed = T)[[1]][3], method_dict)
  for(j in 1:8){
    many_aucs[[j]][[i]] <- many_raw[[j]][many_raw[[j]]$methods_name == best_method,]
    many_aucs[[j]][[i]]$disease <- convert_names(all_authors[i], author_dict)
  }
}

for(j in 1:8){
  many_aucs[[j]] <- do.call("rbind", many_aucs[[j]])
  many_aucs[[j]]$disease <- factor(as.character(many_aucs[[j]]$disease), levels = rev(sort(use_disease)))
}
#need to also read in the new many_auc
#


sex_auc_df <- rbind(sex_auc_df, do.call("rbind", temp_list))
sex_auc_df$disease <- factor(as.character(sex_auc_df$disease), levels = rev(sort(use_disease)))

ethnic_auc_df <- ethnic_auc_df[ethnic_auc_df$best == 1,]
ethnic_auc_df$disease <- factor(as.character(ethnic_auc_df$disease), levels = rev(sort(use_disease)))

sibs_df <- sibs_df[sibs_df$best == 1,]
sibs_df$disease <- factor(as.character(sibs_df$disease), levels = rev(sort(use_disease)))

age_auc_df <- age_auc_df[age_auc_df$best == 1,]
age_auc_df$disease <- factor(as.character(age_auc_df$disease), levels = rev(sort(use_disease)))



#Plotting 
the_plot <- ggplot(score_pheno_def_df, aes(diff, disease, color = variable)) + geom_vline(xintercept = 0) + geom_point() +
  geom_hline(yintercept = 1:length(unique(score_pheno_def_df$disease)) + 0.5, color = "grey80") +
  labs(x = "Diff. From Primary Pheno. AUC", color = "Pheno. Method:", y = "") +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size = 10),
        legend.spacing.x = unit(0, 'cm'),
        legend.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) 
plot(the_plot)
ggsave("meta_plots/meta/for_fig.score.png", the_plot, height = 6, width = 5.5)




the_plot <- ggplot(sex_auc_df, aes(male_auc, disease, color = sex)) + geom_point() +
  geom_hline(yintercept = 1:length(unique(sex_auc_df$disease)) + 0.5, color = "grey80") +
  labs(x = "AUC", y = "", color = "Sex:") +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size = 10),
        legend.spacing.x = unit(0, 'cm'),
        legend.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) 
plot(the_plot)
ggsave("meta_plots/meta/for_fig.sex.png", the_plot, height = 6, width = 5.5)



the_plot <- ggplot(ethnic_auc_df, aes(diff, disease, color = ethnic)) + geom_vline(xintercept = 0) + geom_point() +
  geom_hline(yintercept = 1:length(unique(ethnic_auc_df$disease)) + 0.5, color = "grey80") +
  labs(x = "Diff. from UK PRS", y = "", color = "Population:") +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size = 10),
        legend.spacing.x = unit(0, 'cm'),
        legend.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) 
plot(the_plot)
ggsave("meta_plots/meta/for_fig.ethnic.png", the_plot, height = 6, width = 5.5)



the_plot <- ggplot(sibs_df, aes(conc, disease, color = type)) + geom_point() +
  geom_hline(yintercept = 1:length(unique(sibs_df$disease)) + 0.5, color = "grey80") +
  labs(x = "Concordance", y = "", color = "Relationship:") +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size = 10),
        legend.spacing.x = unit(0, 'cm'),
        legend.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) 
plot(the_plot)
ggsave("meta_plots/meta/for_fig.sib.png", the_plot, height = 6, width = 5.5)




the_plot <- ggplot(age_auc_df, aes(old_auc, disease, color = sex)) + geom_point() +
  geom_hline(yintercept = 1:length(unique(age_auc_df$disease)) + 0.5, color = "grey80") +
  labs(x = "AUC", y = "", color = "Age Group:") +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size = 10),
        legend.spacing.x = unit(0, 'cm'),
        legend.title=element_text(size=12)) +
  theme(plot.margin=unit(c(1,4,1,1),"cm")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) 
plot(the_plot)
ggsave("meta_plots/meta/for_fig.age.png", the_plot, height = 6, width = 5.5)



vars_examined <- c("Time at\nCurrent\nAddress", "Income", "Number in\nHouse", "Year of\nEducation",
                   "Census:\nMedian Age", "Census:\nUnemployment", "Census:\nHealth", "Census:\nPop. Den.")
vars_simple <- c("address", "income", "house_size", "education", "cen_age", "cen_employ", "cen_health", "cen_popdens")
var_splits <- list(c("1-20 years", "20+ years"), c("< 40,000", "> 40,000"), c("1-2 person", "3+ persons"),
                   c("1-19 years", "20+ years"), c("> 42 years", "< 42 years"), c("> 38", "< 38"),
                   c("> 719", "< 719"), c("> 32 persons/hectare", "< 32 persons/hectare"))
for(j in 1:8){
  the_plot <- ggplot(many_aucs[[j]], aes(upper_auc, disease, color = var)) + geom_point() +
    geom_hline(yintercept = 1:length(unique(many_aucs[[j]]$disease)) + 0.5, color = "grey80") +
    labs(x = "AUC", y = "", color = paste0(vars_examined[j], ":")) +
    theme(legend.position = "bottom") +
    theme(legend.text=element_text(size = 10),
          legend.spacing.x = unit(0, 'cm'),
          legend.title=element_text(size=12)) +
    theme(plot.margin=unit(c(1,4,1,1),"cm")) +
    scale_color_discrete(labels = rev(var_splits[[j]])) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) 
  plot(the_plot)
  ggsave(paste0("meta_plots/meta/for_fig.", vars_simple[j], ".png"), the_plot, height = 6, width = 5.5)
  
}


########################################################################################################################
#                                        AUC TEST
#########################################################################################################################

auc_files <- list.files("nonpred_results/per_score/", "auctest")
auc_files <- auc_files[-length(auc_files)]
just_best <- matrix(NA, nrow = length(auc_files), ncol = 2)

for(i in 1:nrow(just_best)){
    curr_file <- readRDS(paste0("nonpred_results/per_score/", auc_files[i]))
    best_method <- read.table(paste0("tune_results/", strsplit(auc_files[i], ".", fixed = T)[[1]][1], ".methods.ss"))
    best_all <- read.table(paste0("tune_results/", strsplit(auc_files[i], ".", fixed = T)[[1]][1], ".best.ss"))
    if(!is.null(curr_file$auc_p$sex)){
      just_best[i, 1] <- curr_file$auc_p$sex[curr_file[["extra"]]$score_names == best_all[3,2]]
    }
    just_best[i, 2] <- curr_file$auc_p$age[curr_file[["extra"]]$score_names == best_all[3,2]]
}

just_best <- data.frame(disease = author_dict[-c(4,23,26),2], just_best, stringsAsFactors = F)
just_best[,2] <- signif(just_best[,2], 3)
just_best[,3] <- signif(just_best[,3], 3)
write.table(just_best, "supp_tables/auc_test_age_sex.txt", row.names = F, col.names = T, quote = F, sep = "\t")

sum(just_best[!is.na(just_best[,2]),2] < 0.05)
sum(just_best[,3] < 0.05)