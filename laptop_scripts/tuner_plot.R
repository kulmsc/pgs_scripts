library(stringr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

author <- "Christophersen"

res <- readRDS(paste0("tune_results/", author, "_res.RDS"))

split_score_names <- str_split(res[["score_names"]], fixed("."), simplify = T)
simple_name <- paste0(split_score_names[,3], "-", split_score_names[,2])
method_name <- split_score_names[,3]


do_plot_series <- function(subset_inds, draw_plot_option = "none", save_option = FALSE, best_stat = "auc"){
  
   if(length(unique(method_name[subset_inds])) == 1){
    plot_ext <- unique(method_name[subset_inds])
   } else {
    plot_ext <- "best"
   }
  
  ##############################################################################
  #Average all of the results
  mean_conc <- Reduce("+", res[["conc"]])/length(res[["conc"]])
  
  mean_survfit <- matrix(0, nrow = length(res[["survfit"]][[1]]), ncol = 6)
  for(i in 1:length(res[["survfit"]])){
    for(j in 1:length(res[["score_names"]])){
      mean_survfit[j,] <- mean_survfit[j,] + as.numeric(tail(res[["survfit"]][[i]][[j]],1)[-1])
    }
  }
  mean_survfit <- mean_survfit/length(res[["survfit"]])
  
  mean_auc <- Reduce("+", res[["auc"]])/length(res[["auc"]])
  
  mean_or <- Reduce("+", res[["or"]])/length(res[["or"]])
  
  
  ###########################################################################
  #Average for the base items
  base_conc <- rep(0, 3)
  for(i in 1:length(res[["base"]][[1]])){
    base_conc <- base_conc + as.numeric(c(res[["base"]][[1]][[i]]$concordance,
                                          res[["base"]][[1]][[i]]$std.err,
                                          res[["base"]][[1]][[i]]$std.err * 1.96))
  }
  base_conc <- base_conc/length(res[["base"]][[1]])
  base_conc_mat <- data.frame("base", "base", t(base_conc), stringsAsFactors = F)
  colnames(base_conc_mat) <- c("score_name", "method", "conc", "se", "ci")
  
  base_survfit <- matrix(0, nrow = 1, ncol = 6)
  for(i in 1:length(res[["base"]][[2]])){
    base_survfit <- base_survfit + as.numeric(tail(res[["base"]][[2]][[i]],1)[-1])
  }
  base_survfit <- base_survfit/9
  
  base_auc <- 0
  for(i in 1:length(res[["base"]][[3]])){
    base_auc <- res[["base"]][[3]][[i]] + base_auc
  }
  base_auc <- base_auc/length(res[["base"]][[3]])
  
  base_or <- 0
  for(i in 1:length(res[["base"]][[4]])){
    base_or <- res[["base"]][[4]][[i]] + base_or
  }
  base_or <- base_or/length(res[["base"]][[4]])

  #################################################################
   #CONC PLOT
   conc <- data.frame(simple_name, method_name, mean_conc)
   conc <- conc[subset_inds,]
   colnames(conc) <- c("score_name", "method", "conc", "se")
   conc$ci <- conc$se * 1.96
   conc$score_name <- factor(conc$score_name, levels = conc$score_name[order(conc$conc)])
   the_conc_plot <- ggplot(conc, aes(conc, score_name)) + geom_point() + 
    geom_errorbarh(aes(xmin = conc-ci, xmax = conc+ci, height = 0.4)) +
    labs(x = "Concordance", y = "Score Name") +
    geom_vline(xintercept = base_conc[1], linetype = "dashed", alpha = 0.3)
  
 
 
   #AUC PLOT
   auc <- data.frame(simple_name, method_name, mean_auc)
   auc <- auc[subset_inds,]
   colnames(auc) <- c("score_name", "method", "ci_lo", "auc", "ci_hi")
   auc$score_name <- factor(auc$score_name, levels = auc$score_name[order(auc$auc)])
   the_auc_plot <- ggplot(auc, aes(auc, score_name)) + geom_point() + 
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, height = 0.4)) +
    labs(x = "AUC", y = "Score Name") +
    geom_vline(xintercept = base_auc[2], linetype = "dashed", alpha = 0.3)
 
 
 
   #EVENT RATE
   km_df <- mean_survfit[subset_inds,]
   km_df <- data.frame(mean_vals = as.numeric(t(km_df[,1:3])), sd_vals = as.numeric(t(km_df[,4:6])))
   km_df$group <- as.factor(rep(1:3, nrow(km_df)/3))
   km_df$simple_name <- rep(simple_name[subset_inds], each = 3)
   km_df$method_name <- rep(method_name[subset_inds], each = 3)
   km_df$simple_name <- factor(km_df$simple_name, levels = km_df$simple_name[km_df$group == 3][
     order(km_df$mean_vals[km_df$group == 3])])

    the_km_plot <- ggplot(km_df, aes(simple_name, mean_vals, color = group,
                                     ymin = mean_vals - sd_vals, ymax = mean_vals + sd_vals)) +
     geom_point(position = position_dodge(width = 0.1)) +
     geom_errorbar(position = position_dodge(width = 0.1), width = 0.2) +
     labs(y = "Cumulative Hazard", x = "Score Name", color = "Risk\nGroup") +
     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      geom_hline(yintercept = base_survfit[1], linetype = "dashed", alpha = 0.3) +
      geom_hline(yintercept = base_survfit[2], linetype = "solid", alpha = 0.3) +
      geom_hline(yintercept = base_survfit[3], linetype = "dashed", alpha = 0.3)
 
 
 
   #ODDS RATIO
   odds_df <- data.frame(mean_or, simple_name, method_name)
   odds_df <- odds_df[subset_inds,]
   colnames(odds_df) <- c("lo", "or", "hi", "simple_name", "method_name")
   odds_df$simple_name <- factor(simple_name[subset_inds], simple_name[subset_inds][order(odds_df$or)])
   the_or_plot <- ggplot(odds_df, aes(or, simple_name)) + geom_point() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.5) +
    labs(x = "Odds Ratio", y = "Score Name") +
     geom_vline(xintercept = base_or[2], linetype = "dashed", alpha = 0.3)

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
}


best_names <- rep("", length(unique(method_name)))
for(i in 1:length(unique(method_name))){
 method_inds <- which(method_name %in% unique(method_name)[i])
 best_names[i] <- do_plot_series(method_inds)
}


final_best_name <- do_plot_series(which(simple_name %in% best_names), "all", TRUE)

