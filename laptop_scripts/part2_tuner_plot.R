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


