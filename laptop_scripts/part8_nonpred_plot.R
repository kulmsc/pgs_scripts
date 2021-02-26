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
  
