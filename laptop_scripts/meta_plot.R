library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(viridis)
theme_set(theme_cowplot())


meta_res <- readRDS("nonpred_results/meta_results.RDS")
meta_convert <- read.table("nonpred_results/meta_name_convert.txt", stringsAsFactors = F, sep = "\t")
method_convert <- read.table("nonpred_results/method_convert.txt", stringsAsFactors = F)
author_convert <- read.table("nonpred_results/author_convert.txt", stringsAsFactors = F, sep = "\t")

authors <- meta_res[["authors"]]
methods <- unique(meta_res[["method_auc"]]$method)

#For making nice disease names, change the old nmaes with underscores into something without
get_new_labels <- function(old_labels, decoder){
  old_labels <- as.character(old_labels)
  new_labels <- rep("", length(old_labels))

  for(i in 1:length(old_labels)){
    if(old_labels[i] != ""){
      new_labels[i] <- decoder[decoder[,1] == old_labels[i], 2]
    }
  }
  return(new_labels)
}


##########################################################################

#Read in all of the meta-statistics and create a few more stats
meta_vals <- meta_res[["best_auc"]]
meta_vals <- meta_vals[,c("uk_control", "uk_case", "uk_snps", "sampe_size", "cases", "controls", "snps", "h2")]
meta_vals$uk_ratio <- meta_vals$uk_case/(meta_vals$uk_control + meta_vals$uk_case)
meta_vals$ss_ratio <- meta_vals$cases/(meta_vals$cases + meta_vals$controls)


#Calculcate pair-wise correlations
corr_vals <- as.data.frame(matrix(0, nrow = ncol(meta_vals), ncol = ncol(meta_vals)))
colnames(corr_vals) <- colnames(meta_vals)
for(i in 1:ncol(meta_vals)){
  for(j in 1:ncol(meta_vals)){
    corr_vals[i,j] <- cor(meta_vals[,i], meta_vals[,j])
  }
}

#Plot the heat map of correlations
corr_vals <- melt(corr_vals)
corr_vals$x <- as.factor(rep(1:ncol(meta_vals), ncol(meta_vals)))
corr_vals$y <- as.factor(rep(1:ncol(meta_vals), each = ncol(meta_vals)))

heatmap_plot <- ggplot(corr_vals, aes(x, y, fill = value)) + geom_raster() +
  scale_x_discrete(labels = get_new_labels(colnames(meta_vals), meta_convert)) +
  scale_y_discrete(labels = get_new_labels(colnames(meta_vals), meta_convert)) +
  labs(x = "", y = "", fill = "Corr.") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_viridis()
plot(heatmap_plot)

####################################################################################


method_coef <- meta_res[["method_coef"]]
method_names <- c("sampe_size", "cases", "cases/controls", "snps", "h2", "uk_snps", "uk_case/uk_control")

#Extract the p-value of the meta-stat feature in each of the meta-regressions
p_vals <- as.data.frame(matrix(0, nrow = length(method_names), ncol = length(method_coef)))
for(i in 1:nrow(p_vals)){
  for(j in 1:ncol(p_vals)){
    p_vals[i,j] <- method_coef[[j]][i,4]
  }
}
colnames(p_vals) <- methods
p_vals$coef <- method_names
p_vals <- melt(p_vals, id.vars = "coef")
p_vals$x <- as.factor(rep(1:length(method_names), length(method_coef)))
p_vals$y <- as.factor(sort(rep(1:length(method_coef), length(method_names))))

#Create a heatmap over each method and meta-stat of the p-values
corr_plot <- ggplot(p_vals, aes(x, y, fill = value)) + geom_raster() +
  scale_x_discrete(labels = get_new_labels(method_names, meta_convert)) +
  scale_y_discrete(labels = get_new_labels(methods, method_convert)) +
  labs(x = "", y = "", fill = "P-Value") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_viridis()
plot(corr_plot)

#Additionaly create just one plot per method for the beta and std. error
all_method_corr_plot <- list()
for(i in 1:length(method_coef)){
  plot_df <- as.data.frame(method_coef[[i]])
  colnames(plot_df) <- c("beta", "se", "z", "p")
  plot_df$stats <- factor(method_names, levels = method_names[order(plot_df$beta)])
  
  method_corr_plot <- ggplot(plot_df, aes(beta, stats)) + geom_point() +
    geom_errorbarh(aes(xmin = beta - se, xmax = beta + se), height = 0.2) +
    labs(x = "Beta", y = "", caption = get_new_labels(methods[i], method_convert)) +
    scale_y_discrete(labels = get_new_labels(levels(plot_df$stats), meta_convert))
  plot(method_corr_plot)
  all_method_corr_plot[[i]] <- method_corr_plot
}

#######################################################################

#Make the same plot as above, but now the regression is over just the best methods for each phenotype
plot_df <- as.data.frame(meta_res[["abs_coef"]])
colnames(plot_df) <- c("beta", "se", "z", "p")
plot_df$stats <- method_names
plot_df$stats <- factor(method_names, levels = method_names[order(plot_df$beta)])

best_method_corr_plot <- ggplot(plot_df, aes(beta, stats)) + geom_point() +
  geom_errorbarh(aes(xmin = beta - se, xmax = beta + se), height = 0.2) +
  labs(x = "Beta", y = "") +
  scale_y_discrete(labels = get_new_labels(levels(plot_df$stats), meta_convert))
plot(best_method_corr_plot)

#######################################################################

#Plot the difference in training and testing AUC
ttcomp <- as.data.frame(meta_res[["ttcomp"]])
colnames(ttcomp) <- c("diff_conc", "diff_survfit", "diff_auc", "diff_or")
ttcomp$author <- authors

tt_plot <- ggplot(ttcomp, aes(diff_auc, authors)) + geom_point() +
  geom_vline(aes(xintercept = 0)) +
  labs(x = "Difference in AUC", y = "Authors") +
  scale_y_discrete(labels = get_new_labels(ttcomp$author, author_convert))
plot(tt_plot)

###########################################################################

for(i in 1:length(method_coef)){
  ggsave(paste0("meta_plots/", methods[i], "_beta_plot.png"), 
         all_method_corr_plot[[i]], "png", width = 5, height = 5)
}

ggsave("meta_plots/best_method_corr_plot.png", best_method_corr_plot, "png", width = 5, height = 5)
ggsave("meta_plots/pval_plot.png", corr_plot, "png", width = 5, height = 5)
ggsave("meta_plots/heatmap_plot.png", heatmap_plot, "png", width = 5, height = 5)
ggsave("meta_plots/tt_plot.png", tt_plot, "png", width = 5, height = 5)
