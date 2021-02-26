#################################################################
##################### ETHNICITY ###############################
#################################################################



ethnic_res <- nonpred_res[["new_ethnic"]]
all_dist <- list()
all_param <- list()
i <- 1
k <- 1
for(ethnic in c("brit_stats", "euro_stats", "african_stats", "asian_stats")){
  for(j in 1:nrow(ethnic_res[[ethnic]])){
    all_dist[[i]] <- data.frame(score = rnorm(1000, mean = ethnic_res[[ethnic]][j,1], sd = ethnic_res[[ethnic]][j,2]))
    all_dist[[i]]$ethnic <- ethnic
    all_dist[[i]]$name <- ethnic_res[["names_to_keep"]][j]

    i <- i + 1
  }
  all_param[[k]] <- data.frame(ethnic_res[[ethnic]][,1:2])
  colnames(all_param[[k]]) <- c("mean", "sd")
  all_param[[k]]$ethnic <- ethnic
  all_param[[k]]$name <- ethnic_res[["names_to_keep"]][1:nrow(all_param[[k]])]
  k <- k + 1
}

all_dist <- do.call("rbind", all_dist) #can create normal distribution looking things from here
all_param <- do.call("rbind", all_param)
all_param$method <- convert_names(str_split(all_param$name, fixed("."), simplify = T)[,3])
for_convert <- cbind(c("brit_stats", "euro_stats", "african_stats", "asian_stats"), c("British", "European", "African", "Asian"))
all_param$ethnic <- convert_names(all_param$ethnic, for_convert)
all_param$method <- factor(all_param$method, levels =
                             all_param$method[all_param$ethnic == "African"][order(all_param$mean[all_param$ethnic == "African"])])

the_plot <- ggplot(all_param, aes(method, mean, color = ethnic)) + geom_point(position=position_dodge(width = 0.9)) + coord_flip() +
  geom_errorbar(position=position_dodge(width = 0.9), aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
  labs(y = "Score", x = "Method", color = "Ethnicity") +
  geom_vline(xintercept = 1:length(unique(all_param$method)) + 0.5, color = "grey80")
plot(the_plot)
ggsave(paste0("meta_plots/", tolower(author), ".nonpred.ethnic.1.png"),
       the_plot, "png", height=6, width=7)
saveRDS(all_param, paste0("nonpred_results/derived_from_per_score/", author, ".ethnic.plot.RDS"))


ethnic_stat <- ethnic_res[["stat_test"]]
keep1 <- ethnic_stat$comp1
keep2 <- ethnic_stat$comp2
colnames(ethnic_stat) <- c(ethnic_res[["names_to_keep"]][-length(ethnic_res[["names_to_keep"]])], "eth1", "eth2")
saveRDS(ethnic_stat, paste0("nonpred_results/derived_from_per_score/", author, ".ethnic.tstat.RDS"))

ethnic_stat <- as.data.frame(ethnic_res[["pval_test"]])
ethnic_stat$comp1 <- keep1
ethnic_stat$comp2 <- keep2
colnames(ethnic_stat) <- c(ethnic_res[["names_to_keep"]][-length(ethnic_res[["names_to_keep"]])], "eth1", "eth2")
saveRDS(ethnic_stat, paste0("nonpred_results/derived_from_per_score/", author, ".ethnic.pval.RDS"))


