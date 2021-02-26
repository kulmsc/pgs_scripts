#################################################################
##################### SIBLING ###############################
#################################################################

sibs_res <- as.data.frame(nonpred_res[["sibs"]][[1]])

sibs_res$method <- convert_names(str_split(rownames(sibs_res), fixed("."), simplify = T)[,3])

sibs_res_sib <- sibs_res[,c(1:3,7)]
colnames(sibs_res_sib) <- c("lo", "conc", "hi", "method")
sibs_res_sib$type <- "Sibling"

sibs_res_non <- sibs_res[,c(4:7)]
colnames(sibs_res_non) <- c("lo", "conc", "hi", "method")
sibs_res_non$type <- "Non-Sib"

sibs_res <- rbind(sibs_res_sib, sibs_res_non)
mean_vals <- unlist(lapply(unique(sibs_res$method), function(x) mean(sibs_res$conc[sibs_res$method == x])))
sibs_res$method <- factor(sibs_res$method, levels = unique(sibs_res$method)[order(mean_vals)])

the_plot <- ggplot(sibs_res, aes(conc, method, color = type)) + geom_boxplot() +
  labs(x = "Concordance", y = "Method", color = "Comparison") +
  geom_hline(yintercept = 1:length(unique(sibs_res$method)) + 0.5, color = "grey80")
plot(the_plot)

ggsave(paste0("meta_plots/", tolower(author), ".nonpred.sibs.1.png"),
       the_plot, "png", height=6, width=7)
saveRDS(sibs_res, paste0("nonpred_results/derived_from_per_score/", author, ".sibs_res.plot.RDS"))


