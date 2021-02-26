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
