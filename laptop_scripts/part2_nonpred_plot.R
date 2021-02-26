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

