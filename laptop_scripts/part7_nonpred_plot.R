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
