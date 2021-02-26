library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(viridis)
library(pROC)
theme_set(theme_cowplot())

author <- "Bentham"
all_authors <- unique(str_split(list.files("nonpred_results/per_score/"), fixed("."), simplify = T)[,1])
all_authors <- all_authors[-which(all_authors == "xie")]

for(author in all_authors){
  print(author)

nonpred_res <- readRDS(paste0("nonpred_results/per_score/", tolower(author), ".res.many.RDS"))
score_names <- nonpred_res[["extra"]][[1]]
method_names <- str_split(score_names, fixed("."), simplify = T)[,3]
best_in_method <- read.table(paste0("tune_results/", tolower(author), ".methods.ss"))

convert_names <- function(x, the_dict = method_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
disease_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = "\t")

neat_method_names <- convert_names(method_names)
neat_score_names <- paste0(neat_method_names, "-", str_split(score_names, fixed("."), simplify = T)[,2])

#################################################################
##################### SCORE SIZES ###############################
#################################################################

#score sizes
score_size_res <- nonpred_res[["score_sizes"]]

all_len <- score_size_res[["all_len"]]
plot_df <- data.frame(neat_score_names, neat_method_names, all_len)
plot_df$neat_score_names <- factor(plot_df$neat_score_names, levels = plot_df$neat_score_names[order(plot_df$all_len)])

the_plot <- ggplot(plot_df, aes(all_len, neat_score_names, color = neat_method_names)) + geom_point() +
  labs(x = "Number SNPs in Score", y = "Score Names", color = "Method") +
  theme(axis.text.y = element_blank())
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.scoresize.1.png"),
       the_plot, "png", height=6, width=7)

mean_vals <- unlist(lapply(as.character(unique(plot_df$neat_method_names)), function(x) mean(plot_df$all_len[plot_df$neat_method_names == x])))
plot_df$neat_method_names <- factor(plot_df$neat_method_names, levels = as.character(unique(plot_df$neat_method_names))[order(mean_vals)])
the_plot <- ggplot(plot_df, aes(log10(all_len), neat_method_names)) + geom_boxplot() + geom_point() +
  labs(x = "log10(Number SNPs in Score)", y = "Method") +
  geom_hline(yintercept = 1:length(unique(plot_df$neat_method_names)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.scoresize.2.png"),
       the_plot, "png", height=6, width=7)


saveRDS(plot_df, paste0("nonpred_results/derived_from_per_score/", author, ".score.len.RDS"))
###

#Equal Splits Mean
score_equal_splits_mean <- data.frame(score_size_res[["all_equal_splits_mean"]])
score_equal_splits_mean <- score_equal_splits_mean/rowSums(score_equal_splits_mean)
colnames(score_equal_splits_mean) <- c("Q1", "Q2", "Q3", "Q4")
score_equal_splits_mean$method <- neat_method_names
score_equal_splits_mean$score <- neat_score_names
save_labels <- score_equal_splits_mean$score[order(score_equal_splits_mean$Q1)]

score_equal_splits_mean <- melt(score_equal_splits_mean, id.vars = c("method", "score"))
score_equal_splits_mean$score <- factor(score_equal_splits_mean$score, levels = save_labels)
the_plot <- ggplot(score_equal_splits_mean, aes(value, score, shape = variable, color = method)) + geom_point() +
  labs(x = "Proportion of Mean Effect", y = "Score", color = "Quartile of\nAbs. Effect") +
  theme(axis.text.y = element_blank())
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.scoremean.1.png"),
       the_plot, "png", height=6, width=7)
saveRDS(score_equal_splits_mean, paste0("nonpred_results/derived_from_per_score/", author, ".score.mean_quant.RDS"))


the_plot <- ggplot(score_equal_splits_mean, aes(value, method, color = variable)) + geom_boxplot() +
  labs(x = "Proportion of Mean Effect", y = "Method", color = "Quartile") +
  geom_hline(yintercept = 1:length(unique(score_equal_splits_mean$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.scoremean.2.png"),
       the_plot, "png", height=6, width=7)

###

#Equal Splits Len
score_equal_splits_len <- data.frame(score_size_res[["all_equal_splits_len"]])
score_equal_splits_len <- score_equal_splits_len/rowSums(score_equal_splits_len)
colnames(score_equal_splits_len) <- c("Q1", "Q2", "Q3", "Q4")
score_equal_splits_len$method <- neat_method_names
score_equal_splits_len$score <- neat_score_names
save_labels <- score_equal_splits_len$score[order(score_equal_splits_len$Q1)]

score_equal_splits_len <- melt(score_equal_splits_len, id.vars = c("method", "score"))
score_equal_splits_len$score <- factor(score_equal_splits_len$score, levels = save_labels)
the_plot <- ggplot(score_equal_splits_len, aes(value, score, shape = variable, color = method)) + geom_point() +
  labs(x = "Proportion of Length", y = "Score", color = "Quartile of\nAbs. Effect") +
  theme(axis.text.y = element_blank())
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.scoresplit.1.png"),
       the_plot, "png", height=6, width=7)
saveRDS(score_equal_splits_len, paste0("nonpred_results/derived_from_per_score/", author, ".score.len_quant.RDS"))

the_plot <- ggplot(score_equal_splits_len, aes(value, method, color = variable)) + geom_boxplot() +
  labs(x = "Proportion of Length", y = "Method", color = "Quartile") +
  geom_hline(yintercept = 1:length(unique(score_equal_splits_len$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.scoresplit.2.png"),
       the_plot, "png", height=6, width=7)

###

#Quantiles
score_quantiles <- data.frame(score_size_res[["all_quantiles"]])
score_quantiles <- score_quantiles/rowSums(score_quantiles)
colnames(score_quantiles) <- paste0("Q.", 1:10)
score_quantiles$method <- neat_method_names
score_quantiles$score <- neat_score_names
score_quantiles$names <- score_names
save_labels <- score_quantiles$names[order(score_quantiles$Q.1)]

score_quantiles <- melt(score_quantiles, id.vars = c("method", "score", "names"))
score_quantiles$names <- factor(score_quantiles$names, levels = save_labels)
score_quantiles$Q <- as.numeric(str_split(score_quantiles$variable, fixed("."), simplify = T)[,2])
the_plot <- ggplot(score_quantiles, aes(value, names, color = Q)) + geom_point() +
  scale_y_discrete(labels = neat_method_names) +
  labs(x = "Proportion of Abs. Effect", y = "Method", color = "Quantile") +
  theme(axis.text=element_text(size=8)) +
  scale_color_viridis()
plot(the_plot)


#################################################################
##################### PHENO DEFS ###############################
#################################################################

pheno_def_res <- nonpred_res[["pheno_defs"]]

#Cases
pheno_def_size <- as.data.frame(pheno_def_res[["total_phens"]])
colnames(pheno_def_size) <- "size"
pheno_def_size$method <- c("ICD", "Self\nReported", "ICD or\nSelf Rep.", "Any", "Double\nReported")
pheno_def_size$method <- factor(pheno_def_size$method, levels = pheno_def_size$method[order(pheno_def_size$size)])
the_plot <- ggplot(pheno_def_size, aes(size, method)) + geom_point() +
  labs(x = "Cases", y = "Phenotyping Method")
saveRDS(pheno_def_size, paste0("nonpred_results/derived_from_per_score/", author, ".pheno_def.len.RDS"))
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.phenodefs.1.png"),
       the_plot, "png", height=6, width=7)

#AUC
pheno_def_auc <- as.data.frame(pheno_def_res[["all_phen_auc"]])
pheno_def_auc <- pheno_def_auc[,seq(2,14,3)]
colnames(pheno_def_auc) <- c("ICD", "Self\nReported", "ICD or\nSelf Rep.", "Any", "Double\nReported")
pheno_def_auc$score <- neat_score_names
pheno_def_auc$method <- neat_method_names

# keep_score_names <- rep("", length(unique(neat_method_names)))
# for(m in unique(neat_method_names)){
#   keep_score_names[which(unique(neat_method_names) == m)] <- pheno_def_auc$score[pheno_def_auc$method == m][
#     which.max(pheno_def_auc$`ICD or\nSelf Rep.`[pheno_def_auc$method == m])]
# }
temp_save <- data.frame(pheno_def_auc)
pheno_def_auc <- melt(pheno_def_auc, id.vars = c("score", "method"))

mean_vals <- unlist(lapply(unique(pheno_def_auc$method),
                           function(x) mean(pheno_def_auc$value[pheno_def_auc$method == x & pheno_def_auc$variable == "ICD or\nSelf Rep."])))
pheno_def_auc$method <- factor(pheno_def_auc$method, levels = unique(pheno_def_auc$method)[order(mean_vals)])
the_plot <- ggplot(pheno_def_auc, aes(value, method, color = variable)) + geom_boxplot() +
  labs(x = "AUC", y = "Method", color = "Phenotyping\nMethod") +
  geom_hline(yintercept = 1:length(unique(pheno_def_auc$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.phenodefs.2.png"),
       the_plot, "png", height=6, width=7)



pheno_def_auc <- temp_save[score_names %in% best_in_method[,1],]
pheno_def_auc <- melt(pheno_def_auc, id.vars = c("score", "method"))
the_plot <- ggplot(pheno_def_auc, aes(variable, method, fill = value)) + geom_raster() + scale_fill_viridis() +
  labs(x = "Phenotyping Method", y = "Method", fill = "AUC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.phenodefs.3.png"),
       the_plot, "png", height=6, width=7)
saveRDS(pheno_def_auc, paste0("nonpred_results/derived_from_per_score/", author, ".pheno_def.auc.RDS"))

#################################################################
##################### SEX DIFF ###############################
#################################################################

if(!is.null(nonpred_res[["sex_split"]][[1]])){
    sex_res <- as.data.frame(nonpred_res[["sex_split"]][[1]])

    # keep_score_names <- rep("", length(unique(neat_method_names)))
    # for(m in unique(neat_method_names)){
    #   keep_score_names[which(unique(neat_method_names) == m)] <-
    #     neat_score_names[neat_method_names == m][which.max(sex_res[neat_method_names == m,2])]
    # }

    colnames(sex_res) <- c("female_lo", "female_auc", "female_hi", "male_lo", "male_auc", "male_hi")
    sex_res$score_name <- neat_score_names
    sex_res$methods_name <- neat_method_names
    temp_save <- data.frame(sex_res)

    plot_df <- melt(sex_res, id.vars = c("score_name", "methods_name"))
    plot_df <- plot_df[plot_df$variable %in% c("female_auc", "male_auc"),]

    mean_vals <- unlist(lapply(unique(plot_df$methods_name), function(x) mean(plot_df$value[plot_df$methods_name == x])))
    plot_df$methods_name <- factor(plot_df$methods_name, levels = unique(plot_df$methods_name)[order(mean_vals)])
    the_plot <- ggplot(plot_df, aes(value, methods_name, color = variable)) + geom_boxplot() +
      labs(x = "AUC", y = "Methods", color = "Sex") +
      scale_color_discrete(labels = c("Female", "Male")) +
      geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80")
    plot(the_plot)
    ggsave(paste0("meta_plots/", tolower(author), ".nonpred.sexdiff.1.png"),
           the_plot, "png", height=6, width=7)


    plot_df <- temp_save
    plot_df <- sex_res[,c(1:3,7,8)]
    colnames(plot_df) <- c("male_lo", "male_auc", "male_hi", "score_name", "methods_name")
    plot_df <- rbind(plot_df, temp_save[,4:8])
    plot_df$sex <- c(rep("female", nrow(sex_res)), rep("male", nrow(sex_res)))

    plot_df <- plot_df[score_names %in% best_in_method[,1],]
    plot_df$sex <- str_to_title(plot_df$sex)
    plot_df$methods_name <- factor(plot_df$methods_name, levels =
                                     plot_df$methods_name[plot_df$sex == "Male"][order(plot_df$male_auc[plot_df$sex == "Male"])])

    the_plot <- ggplot(plot_df, aes(male_auc, methods_name, color = sex)) +
      geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
      geom_errorbarh(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
                     aes(xmin = male_lo, xmax = male_hi, height = 0)) +
      geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80") +
      labs(x = "AUC", color = "Sex", y = "Methods")
    plot(the_plot)
    ggsave(paste0("meta_plots/", tolower(author), ".nonpred.sexdiff.2.png"),
           the_plot, "png", height=6, width=7)
    saveRDS(plot_df, paste0("nonpred_results/derived_from_per_score/", author, ".sex.auc.RDS"))
}

#################################################################
##################### ETHNICITY ###############################
#################################################################



ethnic_res <- nonpred_res[["new_ethnic"]]
all_dist <- list()
all_param <- list()
i <- 1
k <- 1
for(ethnic in c("brit_stats", "euro_stats", "african_stats", "asian_stats")){
  for(j in 1:nrow(ethnic_res[[ethnic]])){
    all_dist[[i]] <- data.frame(score = rnorm(1000, mean = ethnic_res[[ethnic]][j,1], sd = ethnic_res[[ethnic]][j,2]))
    all_dist[[i]]$ethnic <- ethnic
    all_dist[[i]]$name <- ethnic_res[["names_to_keep"]][j]

    i <- i + 1
  }
  all_param[[k]] <- data.frame(ethnic_res[[ethnic]][,1:2])
  colnames(all_param[[k]]) <- c("mean", "sd")
  all_param[[k]]$ethnic <- ethnic
  all_param[[k]]$name <- ethnic_res[["names_to_keep"]][1:nrow(all_param[[k]])]
  k <- k + 1
}

all_dist <- do.call("rbind", all_dist) #can create normal distribution looking things from here
all_param <- do.call("rbind", all_param)
all_param$method <- convert_names(str_split(all_param$name, fixed("."), simplify = T)[,3])
for_convert <- cbind(c("brit_stats", "euro_stats", "african_stats", "asian_stats"), c("British", "European", "African", "Asian"))
all_param$ethnic <- convert_names(all_param$ethnic, for_convert)
all_param$method <- factor(all_param$method, levels =
                             all_param$method[all_param$ethnic == "African"][order(all_param$mean[all_param$ethnic == "African"])])

the_plot <- ggplot(all_param, aes(method, mean, color = ethnic)) + geom_point(position=position_dodge(width = 0.9)) + coord_flip() +
  geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
  labs(y = "Score", x = "Method", color = "Ethnicity") +
  geom_vline(xintercept = 1:length(unique(all_param$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.ethnic.1.png"),
       the_plot, "png", height=6, width=7)
saveRDS(all_param, paste0("nonpred_results/derived_from_per_score/", author, ".ethnic.plot.RDS"))


ethnic_stat <- ethnic_res[["stat_test"]]
keep1 <- ethnic_stat$comp1
keep2 <- ethnic_stat$comp2
colnames(ethnic_stat) <- c(ethnic_res[["names_to_keep"]][-length(ethnic_res[["names_to_keep"]])], "eth1", "eth2")
saveRDS(ethnic_stat, paste0("nonpred_results/derived_from_per_score/", author, ".ethnic.tstat.RDS"))

ethnic_stat <- as.data.frame(ethnic_res[["pval_test"]])
ethnic_stat$comp1 <- keep1
ethnic_stat$comp2 <- keep2
colnames(ethnic_stat) <- c(ethnic_res[["names_to_keep"]][-length(ethnic_res[["names_to_keep"]])], "eth1", "eth2")
saveRDS(ethnic_stat, paste0("nonpred_results/derived_from_per_score/", author, ".ethnic.pval.RDS"))



#################################################################
##################### SIBLING ###############################
#################################################################

sibs_res <- as.data.frame(nonpred_res[["sibs"]][[1]])

sibs_res$method <- convert_names(str_split(rownames(sibs_res), fixed("."), simplify = T)[,3])

sibs_res_sib <- sibs_res[,c(1:3,7)]
colnames(sibs_res_sib) <- c("lo", "conc", "hi", "method")
sibs_res_sib$type <- "Sibling"

sibs_res_non <- sibs_res[,c(4:7)]
colnames(sibs_res_non) <- c("lo", "conc", "hi", "method")
sibs_res_non$type <- "Non-Sib"

sibs_res <- rbind(sibs_res_sib, sibs_res_non)
mean_vals <- unlist(lapply(unique(sibs_res$method), function(x) mean(sibs_res$conc[sibs_res$method == x])))
sibs_res$method <- factor(sibs_res$method, levels = unique(sibs_res$method)[order(mean_vals)])

the_plot <- ggplot(sibs_res, aes(conc, method, color = type)) + geom_boxplot() +
  labs(x = "Concordance", y = "Method", color = "Comparison") +
  geom_hline(yintercept = 1:length(unique(sibs_res$method)) + 0.5, color = "grey80")
plot(the_plot)

ggsave(paste0("meta_plots/", tolower(author), ".nonpred.sibs.1.png"),
       the_plot, "png", height=6, width=7)
saveRDS(sibs_res, paste0("nonpred_results/derived_from_per_score/", author, ".sibs_res.plot.RDS"))




#################################################################
##################### AGE DIFF ###############################
#################################################################



  age_res <- as.data.frame(nonpred_res[["age_diff"]])

  colnames(age_res) <- c("young_lo", "young_auc", "young_hi", "old_lo", "old_auc", "old_hi")
  age_res$score_name <- neat_score_names
  age_res$methods_name <- neat_method_names
  temp_save <- data.frame(age_res)

  plot_df <- melt(age_res, id.vars = c("score_name", "methods_name"))
  plot_df <- plot_df[plot_df$variable %in% c("young_auc", "old_auc"),]

  mean_vals <- unlist(lapply(unique(plot_df$methods_name), function(x) mean(plot_df$value[plot_df$methods_name == x])))
  plot_df$methods_name <- factor(plot_df$methods_name, levels = unique(plot_df$methods_name)[order(mean_vals)])
  the_plot <- ggplot(plot_df, aes(value, methods_name, color = variable)) + geom_boxplot() +
    labs(x = "AUC", y = "Methods", color = "Age\nGroup") +
    scale_color_discrete(labels = c("Young", "Old")) +
    geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80")
  plot(the_plot)
  ggsave(paste0("meta_plots/", tolower(author), ".nonpred.agediff.1.png"),
         the_plot, "png", height=6, width=7)


  plot_df <- temp_save
  plot_df <- age_res[,c(1:3,7,8)]
  colnames(plot_df) <- c("old_lo", "old_auc", "old_hi", "score_name", "methods_name")
  plot_df <- rbind(plot_df, temp_save[,4:8])
  plot_df$sex <- c(rep("young", nrow(age_res)), rep("old", nrow(age_res)))

  plot_df <- plot_df[score_names %in% best_in_method[,1],]
  plot_df$sex <- str_to_title(plot_df$sex)
  plot_df$methods_name <- factor(plot_df$methods_name, levels =
                                   plot_df$methods_name[plot_df$sex == "Old"][order(plot_df$old_auc[plot_df$sex == "Old"])])
  the_plot <- ggplot(plot_df, aes(old_auc, methods_name, color = sex)) +
    geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
    geom_errorbarh(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
                   aes(xmin = old_lo, xmax = old_hi, height = 0)) +
    geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80") +
    labs(x = "AUC", color = "Age\nGroup", y = "Methods")
  plot(the_plot)
  ggsave(paste0("meta_plots/", tolower(author), ".nonpred.agediff.2.png"),
         the_plot, "png", height=6, width=7)
  saveRDS(plot_df, paste0("nonpred_results/derived_from_per_score/", author, ".agediff.auc.RDS"))
  
  
  #################################################################
  ##################### MANY OTHER AUC ###############################
  #################################################################
  


  vars_examined <- c("Time at\nCurrent\nAddress", "Income", "Number in\nHouse", "Year of\nEducation",
                     "Census:\nMedian Age", "Census:\nUnemployment", "Census:\nHealth", "Census:\nPop. Den.")
  var_splits <- list(c("1-20 years", "20+ years"), c("< 40,000", "> 40,000"), c("1-2 person", "3+ persons"),
                     c("1-19 years", "20+ years"), c("> 42 years", "< 42 years"), c("> 38", "< 38"),
                     c("> 719", "< 719"), c("> 32 persons/hectare", "< 32 persons/hectare"))
  all_plot_df <- list()
  all_pval_df <- list()
  
  for(var_ind in 1:length(vars_examined)){
    
    new_res <- as.data.frame(nonpred_res[["many_auc"]][[var_ind]])
    new_pval <- as.data.frame(nonpred_res[["many_pval"]][[var_ind]])
    
    colnames(new_res) <- c("bottom_lo", "bottom_auc", "bottom_hi", "upper_lo", "upper_auc", "upper_hi")
    new_res$score_name <- neat_score_names
    new_res$methods_name <- neat_method_names
    temp_save <- data.frame(new_res)
    
    plot_df <- melt(new_res, id.vars = c("score_name", "methods_name"))
    plot_df <- plot_df[plot_df$variable %in% c("bottom_auc", "upper_auc"),]
    
    mean_vals <- unlist(lapply(unique(plot_df$methods_name), function(x) mean(plot_df$value[plot_df$methods_name == x])))
    plot_df$methods_name <- factor(plot_df$methods_name, levels = unique(plot_df$methods_name)[order(mean_vals)])
    the_plot <- ggplot(plot_df, aes(value, methods_name, color = variable)) + geom_boxplot() +
      labs(x = "AUC", y = "Methods", color = vars_examined[var_ind]) +
      scale_color_discrete(labels = var_splits[[var_ind]]) +
      geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80")
    plot(the_plot)
    ggsave(paste0("meta_plots/", tolower(author), ".nonpred.many.", var_ind, ".1.png"),
           the_plot, "png", height=6, width=7)
    
    
    plot_df <- temp_save
    plot_df <- new_res[,c(1:3,7,8)]
    colnames(plot_df) <- c("upper_lo", "upper_auc", "upper_hi", "score_name", "methods_name")
    plot_df <- rbind(plot_df, temp_save[,4:8])
    plot_df$var <- c(rep("young", nrow(new_res)), rep("old", nrow(new_res)))
    
    plot_df <- plot_df[score_names %in% best_in_method[,1],]
  
    plot_df$methods_name <- factor(plot_df$methods_name, levels = 
                                     plot_df$methods_name[plot_df$var == "old"][order(plot_df$upper_auc[plot_df$var == "old"])])
    
    the_plot <- ggplot(plot_df, aes(upper_auc, methods_name, color = var)) + 
      geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5)) +
      geom_errorbarh(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5), 
                     aes(xmin = upper_lo, xmax = upper_hi, height = 0)) +
      geom_hline(yintercept = 1:length(unique(plot_df$methods_name)) + 0.5, color = "grey80") +
      labs(x = "AUC", color = vars_examined[var_ind], y = "Methods") +
      scale_color_discrete(labels = rev(var_splits[[var_ind]]))
    plot(the_plot)
    
    ggsave(paste0("meta_plots/", tolower(author), ".nonpred.many.", var_ind, ".2.png"),
           the_plot, "png", height=6, width=7)
    
    new_pval <- new_pval[score_names %in% best_in_method[,1],]
    plot_df$pval <- new_pval
    all_plot_df[[var_ind]] <- plot_df
    
    #new_pval$score_name <- neat_score_names
    #new_pval$methods_name <- neat_method_names
    #new_pval <- new_pval[score_names %in% best_in_method[,1],]
    #all_pval_df[[var_ind]] <- new_pval
    
  }
  
  saveRDS(all_plot_df, paste0("nonpred_results/derived_from_per_score/", author, ".many.auc.RDS"))
  
}
