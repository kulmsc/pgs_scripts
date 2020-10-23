library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(pROC)
theme_set(theme_cowplot())

author <- "Bentham"

nonpred_res <- readRDS(paste0("nonpred_results/", tolower(author), ".res.RDS"))
nonpred_res[[5]] <- readRDS("nonpred_results/extra.RDS")
names(nonpred_res)[5] <- "extra"
score_names <- nonpred_res[["extra"]][[1]]
method_names <- str_split(score_names, fixed("."), simplify = T)[,3]

convert_names <- read.table("nonpred_results/method_convert.txt", stringsAsFactors = F, sep = "\t")

change_names <- function(x){
  for(i in 1:length(convert_names[,1])){
    x[x == convert_names[i,1]] <- convert_names[i,2] 
  }
  return(x)
}

neat_method_names <- change_names(method_names)
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
ggsave(paste0("output_plots/", tolower(author), ".nonpred.scoresize.1.png"),
       the_plot, "png", height=6, width=7)

the_plot <- ggplot(plot_df, aes(log10(all_len), neat_method_names)) + geom_boxplot() + geom_point() +
  labs(x = "log10(Number SNPs in Score)", y = "Method")
plot(the_plot)
ggsave(paste0("output_plots/", tolower(author), ".nonpred.scoresize.2.png"),
       the_plot, "png", height=6, width=7)

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
ggsave(paste0("output_plots/", tolower(author), ".nonpred.scoremean.1.png"),
       the_plot, "png", height=6, width=7)


the_plot <- ggplot(score_equal_splits_mean, aes(value, method, fill = variable)) + geom_boxplot() +
  labs(x = "Proportion of Mean Effect", y = "Method", fill = "Quartile")
plot(the_plot)
ggsave(paste0("output_plots/", tolower(author), ".nonpred.scoremean.2.png"),
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
ggsave(paste0("output_plots/", tolower(author), ".nonpred.scoresplit.1.png"),
       the_plot, "png", height=6, width=7)

the_plot <- ggplot(score_equal_splits_len, aes(value, method, fill = variable)) + geom_boxplot() +
  labs(x = "Proportion of Length", y = "Method", fill = "Quartile")
plot(the_plot)
ggsave(paste0("output_plots/", tolower(author), ".nonpred.scoresplit.2.png"),
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



#################################################################
##################### PHENO DEFS ###############################
#################################################################

pheno_def_res <- nonpred_res[["pheno_defs"]]

#Cases
pheno_def_size <- as.data.frame(pheno_def_res[["total_phens"]])
colnames(pheno_def_size) <- "size"
pheno_def_size$method <- c("ICD", "Self-Rep", "Diag", "Any", "Double")
the_plot <- ggplot(pheno_def_size, aes(size, method)) + geom_point() +
  labs(x = "Cases", y = "Phenotyping Method")
ggsave(paste0("output_plots/", tolower(author), ".nonpred.phenodefs.1.png"),
       the_plot, "png", height=6, width=7)

#AUC
pheno_def_auc <- as.data.frame(pheno_def_res[["all_phen_auc"]])
pheno_def_auc <- pheno_def_auc[,seq(2,14,3)]
colnames(pheno_def_auc) <- c("ICD", "Self-Rep", "Diag", "Any", "Double")
pheno_def_auc$score <- neat_score_names
pheno_def_auc$method <- neat_method_names

keep_score_names <- rep("", length(unique(neat_method_names)))
for(m in unique(neat_method_names)){
  keep_score_names[which(unique(neat_method_names) == m)] <- pheno_def_auc$score[pheno_def_auc$method == m][
    which.max(pheno_def_auc$Diag[pheno_def_auc$method == m])]
}
pheno_def_auc <- melt(pheno_def_auc, id.vars = c("score", "method"))

the_plot <- ggplot(pheno_def_auc, aes(value, method, fill = variable)) + geom_boxplot() +
  labs(x = "AUC", y = "Method", fill = "Phenotyping\nMethod")
plot(the_plot)
ggsave(paste0("output_plots/", tolower(author), ".nonpred.phenodefs.2.png"),
       the_plot, "png", height=6, width=7)

pheno_def_auc <- pheno_def_auc[pheno_def_auc$score %in% keep_score_names,]
the_plot <- ggplot(pheno_def_auc, aes(variable, method, fill = value)) + geom_raster() + scale_fill_viridis() +
  labs(x = "Phenotyping Method", y = "Method", fill = "AUC") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
plot(the_plot)
ggsave(paste0("output_plots/", tolower(author), ".nonpred.phenodefs.3.png"),
       the_plot, "png", height=6, width=7)


#################################################################
##################### ETHNICITY ###############################
#################################################################

ethnic_res <- nonpred_res[["ethnic"]]
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
all_param$method <- change_names(str_split(all_param$name, fixed("."), simplify = T)[,3])

the_plot <- ggplot(all_param, aes(method, mean, color = ethnic)) + geom_point(position=position_dodge(width = 0.9)) + coord_flip() +
  geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
  labs(y = "Score", x = "Method", color = "Ethnicity") 
plot(the_plot)
ggsave(paste0("output_plots/", tolower(author), ".nonpred.ethnic.1.png"),
       the_plot, "png", height=6, width=7)


#################################################################
##################### SCORE SIZES ###############################
#################################################################

sibs_res <- as.data.frame(nonpred_res[["sibs"]][[1]])

sibs_res$method <- change_names(str_split(rownames(sibs_res), fixed("."), simplify = T)[,3])

sibs_res_sib <- sibs_res[,c(1:3,7)]
colnames(sibs_res_sib) <- c("lo", "conc", "hi", "method")
sibs_res_sib$type <- "Sibling"

sibs_res_non <- sibs_res[,c(4:7)]
colnames(sibs_res_non) <- c("lo", "conc", "hi", "method")
sibs_res_non$type <- "Non-Sibg"

sibs_res <- rbind(sibs_res_sib, sibs_res_non)

the_plot <- ggplot(sibs_res, aes(conc, method, fill = type)) + geom_boxplot() +
  labs(x = "Concordance", y = "Method", fill = "Comparison")
plot(the_plot)

ggsave(paste0("output_plots/", tolower(author), ".nonpred.sibs.1.png"),
       the_plot, "png", height=6, width=7)
