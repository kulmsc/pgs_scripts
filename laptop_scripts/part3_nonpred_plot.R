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
