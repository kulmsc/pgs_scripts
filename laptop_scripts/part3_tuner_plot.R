#Get best name among the subsetinds
#################################################################
if(best_stat == "conc"){
  best_name <- as.character(conc$score_name[which.max(conc$conc)])
}else if(best_stat == "auc"){
  best_name <- as.character(auc$score_name[which.max(auc$auc)])
}else if(best_stat == "km"){
  comp_km <- data.frame(name = km_df$simple_name[km_df$group == 1],
                        ratio = km_df$mean_vals[km_df$group == 3]/km_df$mean_vals[km_df$group == 1])
  comp_km$ratio[is.na(comp_km$ratio)] <- 0
  best_name <- as.character(comp_km$name[which.max(comp_km$ratio)])
}else if(best_stat == "or"){
  best_name <- as.character(odds_df$simple_name[which.max(odds_df$or)])
}


#make the plots
#################################################################
if(draw_plot_option == "all"){
  plot(the_auc_plot)
  plot(the_conc_plot)
  plot(the_or_plot)
  plot(the_km_plot)
  
  if(save_option){
    ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".conc.tune.png"),
           the_conc_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".auc.tune.png"),
           the_auc_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".or.tune.png"),
           the_or_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".km.tune.png"),
           the_km_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 7, height = 6)
  }
  
} else if(draw_plot_option == "choose_stat"){
  if(best_stat == "conc"){
    plot(the_conc_plot)
    if(save_option){
      ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".conc.tune.png"),
             the_conc_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    }
    
  }else if(best_stat == "auc"){
    plot(the_auc_plot)
    if(save_option){
      ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".auc.tune.png"),
             the_auc_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    }
    
  }else if(best_stat == "km"){
    plot(the_km_plot)
    if(save_option){
      ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".km.tune.png"),
             the_km_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 7, height = 6)
    }
    
  }else if(best_stat == "or"){
    plot(the_or_plot)
    if(save_option){
      ggsave(paste0("output_plots/", tolower(author), ".", plot_ext, ".or.tune.png"),
             the_or_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), "png", width = 5, height = 5)
    }
    
  }
}

return(best_name)
#}


best_names <- rep("", length(unique(method_name)))
for(i in 1:length(unique(method_name))){
  method_inds <- which(method_name %in% unique(method_name)[i])
  best_names[i] <- do_plot_series(method_inds)
}


final_best_name <- do_plot_series(which(simple_name %in% best_names), "all", TRUE)


