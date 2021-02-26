library(reshape2)
library(ggplot2)
library(cowplot)
library(stringr)
library(viridis)
theme_set(theme_cowplot())


#just need to add method_holder

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
author_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = '\t')

convert_names <- function(x, the_dict = author_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}


#1 - authors
#2 - the stats for best methods for each authors (setup data)
#3 - the stats for the single best methods for each authors (setup data)

meta_results <- readRDS("nonpred_results/meta_results.RDS")

#1 - authors
#2 - methods_auc - all aucs for all diseases and best methods
#3 - best_auc - all aucs for all disease and abs best
#4 - method_coef - regression coefs stratified by best methods
#5 - abs_coef - regression coefs between abs best
#6 - all_coef - regression coefs between all coefs
#7 - method_holder - count the number of times a method appears in the top
#8 - ttcomp

methods_setup <- meta_results[[2]] #nothing to do unless want to plot the AUC vs the meta-stat for each method or author
best_setup <- meta_results[[3]] #can plot the AUC vs meta-stats

temp <- meta_results[["method_auc"]]
temp <- temp[!duplicated(temp$author),]
temp <- temp[,c(6, 12:21, 23)]
temp[,1] <- convert_names(tolower(temp$author), author_dict)
temp <- temp[,-which(colnames(temp) == c("cases", "controls"))]
for(i in 2:ncol(temp)){
  temp[,i] <- signif(temp[,i], 3)
}
write.table(temp, "supp_tables/meta_reg_input.txt", col.names = T, row.names = F, sep = "\t", quote = F)

df <- read.table("supp_tables/med_endpoint_2.txt", stringsAsFactors = F, header = T, sep = "\t", fill = T)
df[,1] <- convert_names(tolower(df$author), author_dict)
write.table(df, "supp_tables/med_endpoint_3.txt", row.names = F, col.names = T, quote = F, sep = "\t")

######################################################################
#                   METHOD REGRESSIONS                            @
######################################################################

reg_coefs <- meta_results[[4]] #list elements is by methods, rows in the element is meta-stat
#c("sampe_size", "cases", "cases/controls", "snps", "h2", "uk_snps", "uk_case/uk_control", "uk_case")



stat_names <- c("uk_snps", "uk_beta_pdiff", "uk_len_pdiff", "sig6_gwas_snps", "sig8_gwas_snps", "len_pdiff",
                "effect_pdiff", "sampe_size", "ancestry", "snps", "h2", "uk_cc_ratio", "gwas_cc_ratio")
plot_stat_names <- c("PRS SNPs", "PRS Dist. by Effect", "PRS Dist. by Length", "GWAS SNPs < 1e-6", "GWAS SNPs < 1e-8",
                     "GWAS Dist. by Length", "GWAS Dist. by Effect", "GWAS Sample Size", "GWAS Ancestry", "GWAS SNPs",
                     "GWAS Heritability", "UKBB CC Ratio", "GWAS CC Ratio")


plot_df <- data.frame(matrix(0, nrow = length(reg_coefs), ncol = 4))
colnames(plot_df) <- c("method", "coef", "se", "p")
poss_methods <- convert_names(unique(methods_setup$method), method_dict)

all_plot_df <- list()
for(stat_ind in 1:length(stat_names)){
  for(i in 1:nrow(plot_df)){
    plot_df[i,1] <- poss_methods[i]
    plot_df[i,2:4] <- reg_coefs[[i]][stat_ind,c(1:2, 4)]
  }
  
  plot_df$method <- factor(plot_df$method, levels = plot_df$method[order(plot_df$coef)])
  
  the_plot <- ggplot(plot_df, aes(coef, method)) + geom_point() +
    geom_errorbar(aes(xmin = coef - se, xmax = coef + se, width = 0)) +
    geom_vline(xintercept = 0) +
    labs(x = "Regression Coefficient", y = "", caption = plot_stat_names[stat_ind]) 
    
  plot(the_plot)
  ggsave(paste0("reg_plots/reg_coef.", tolower(gsub("/", "", gsub(" ",  "_", stat_names)))[stat_ind], ".png"),
         the_plot + theme(plot.margin=unit(c(1,1,1,1), "cm")), height = 6, width = 6)
 
  plot_df$stat_name <- plot_stat_names[stat_ind]
  all_plot_df[[stat_ind]] <- plot_df 
}

plot_df <- do.call("rbind", all_plot_df)
the_plot <- ggplot(plot_df, aes(coef, method)) + geom_point() +
  facet_grid(cols = vars(stat_name), scales = "free") +
  geom_vline(xintercept = 0) +
  labs(x = "Regression Coefficient", y = "") +
  theme(strip.text.x = element_text(size = 8)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 2))
plot(the_plot)
ggsave(paste0("reg_plots/reg_coef.facet.png"), the_plot, height = 6, width = 16)

plot_df[plot_df$p < 0.05,]


mat_beta <- matrix(0, nrow = length(reg_coefs), ncol = length(stat_names))
mat_pval <- matrix(0, nrow = length(reg_coefs), ncol = length(stat_names))
for(i in 1:length(reg_coefs)){
  for(j in 1:length(stat_names)){
    mat_beta[i,j] <- signif(reg_coefs[[i]][j,1], 3)
    mat_pval[i,j] <- signif(reg_coefs[[i]][j,4], 3)
  }
}
mat_beta <- cbind(poss_methods, data.frame(mat_beta))
mat_pval <- cbind(poss_methods, data.frame(mat_pval))
colnames(mat_beta) <- c("Method", plot_stat_names)
colnames(mat_pval) <- c("Method", plot_stat_names)
write.table(mat_beta, "supp_tables/meta_method.beta.txt", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(mat_pval, "supp_tables/meta_method.pval.txt", row.names = F, col.names = T, sep = "\t", quote = F)

write.table(t(mat_beta)[,c(1:8)], "supp_tables/meta_method.beta.1.txt", row.names = T, col.names = F, sep = "\t", quote = F)
write.table(t(mat_beta)[,c(1,9:15)], "supp_tables/meta_method.beta.2.txt", row.names = T, col.names = F, sep = "\t", quote = F)
write.table(t(mat_pval)[,c(1:8)], "supp_tables/meta_method.pval.1.txt", row.names = T, col.names = F, sep = "\t", quote = F)
write.table(t(mat_pval)[,c(1,9:15)], "supp_tables/meta_method.pval.txt.2.txt", row.names = T, col.names = F, sep = "\t", quote = F)


######################################################################
#                   ABS BEST REGRESSIONS                            @
######################################################################

reg_coefs <- data.frame(meta_results[[5]])
colnames(reg_coefs) <- c("coef", "se", "z", "p")
reg_coefs$stat_names <- plot_stat_names
reg_coefs$stat_names <- factor(reg_coefs$stat_names, levels = reg_coefs$stat_names[order(reg_coefs$coef)])
reg_coefs <- reg_coefs[abs(reg_coefs$coef) < 5, ]

the_plot <- ggplot(reg_coefs, aes(coef, stat_names)) + geom_point() +
  geom_errorbar(aes(xmin = coef - se, xmax = coef + se, width = 0)) +
  geom_vline(xintercept = 0) +
  labs(x = "Regression Coefficient", y = "", caption = "Among Best Scores")
plot(the_plot)
ggsave(paste0("reg_plots/reg_coef.best.png"), the_plot, height = 6, width = 7)

for_table <- reg_coefs[,c(5,1,2,4)]
for_table[,2] <- signif(reg_coefs[,2], 3)
for_table[,3] <- signif(reg_coefs[,3], 3)
for_table[,4] <- signif(reg_coefs[,4], 3)
write.table(for_table, "supp_tables/meta_reg.disease.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#to double check
df <- meta_results[["method_auc"]]
best_names <- readRDS("nonpred_results/best_names.RDS")
df <- df[df$name %in% best_names,]

######################################################################
#                   ALL AUC REGRESSIONS                            @
######################################################################

reg_coefs <- data.frame(meta_results[[6]])
colnames(reg_coefs) <- c("coef", "se", "z", "p")
reg_coefs$stat_names <- plot_stat_names
reg_coefs$stat_names <- factor(reg_coefs$stat_names, levels = reg_coefs$stat_names[order(reg_coefs$coef)])
reg_coefs <- reg_coefs[abs(reg_coefs$coef) < 5, ]

the_plot <- ggplot(reg_coefs, aes(coef, stat_names)) + geom_point() +
  geom_errorbar(aes(xmin = coef - se, xmax = coef + se, width = 0)) +
  geom_vline(xintercept = 0) +
  labs(x = "Regression Coefficient", y = "", caption = "Among All Scores")
plot(the_plot)
ggsave(paste0("reg_plots/reg_coef.all.png"), the_plot, height = 6, width = 7)

for_table <- reg_coefs[,c(5,1,2,4)]
for_table[,2] <- signif(reg_coefs[,2], 3)
for_table[,3] <- signif(reg_coefs[,3], 3)
for_table[,4] <- signif(reg_coefs[,4], 3)
write.table(for_table, "supp_tables/meta_reg.all.txt", row.names = F, col.names = T, sep = "\t", quote = F)


######################################################################
#                   METHODS BEST AT EXTREMES                            @
######################################################################
exit()
method_ex <- data.frame(meta_results[[7]])
method_ex$stat <- paste(rep(plot_stat_names[plot_stat_names != "GWAS Ancestry"], each = 2), c("", "Rev."))
method_ex  <- melt(method_ex, id.vars = "stat")
method_ex$variable <- convert_names(method_ex$variable, method_dict)

the_plot <- ggplot(method_ex, aes(variable, stat, fill = value)) + geom_raster() +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "", y = "", fill = "Times In\nTop 10")
plot(the_plot)
ggsave(paste0("reg_plots/method_extreme.png"), the_plot, height = 7, width = 7)

method_ex <- method_ex[grepl("GWAS", method_ex$stat),]
method_ex$stat <- unlist(lapply(strsplit(method_ex$stat, "GWAS"), function(x) x[2]))
method_ex$stat <- trimws(method_ex$stat)
method_ex$stat[grepl("SNP", method_ex$stat)] <- paste0("No. ", method_ex$stat[grepl("SNP", method_ex$stat)])
the_plot <- ggplot(method_ex, aes(variable, stat, fill = value)) + geom_raster() +
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "", y = "", fill = "Times In\nTop 10")
plot(the_plot)
ggsave(paste0("reg_plots/method_extreme_gwas.png"), the_plot, height = 5, width = 7)

to_save <-  data.frame(meta_results[[7]])
to_save$stat <- paste(rep(plot_stat_names[plot_stat_names != "GWAS Ancestry"], each = 2), c("", "Rev."))
colnames(to_save) <- convert_names(colnames(to_save), method_dict)
to_save <- to_save[,c(ncol(to_save), 1:(ncol(to_save)-1))]
write.table(to_save[,c(1,2:7)], "supp_tables/method_extreme.1.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(to_save[,c(1,8:16)], "supp_tables/method_extreme.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(to_save[,c(1,11:16)], "supp_tables/method_extreme.3.txt", row.names = F, col.names = T, quote = F, sep = "\t")

######################################################################
#                  TRAIN/TEST                          @
######################################################################

tt_df <- data.frame(meta_results[[8]])
colnames(tt_df) <- c("diff_conc", "diff_survfit", "diff_auc", "diff_or")
tt_df$disease <- convert_names(meta_results[[1]])

mean(tt_df$diff_auc)
sd(tt_df$diff_auc)/sqrt(nrow(tt_df))

melt_tt <- melt(tt_df, id.vars = "disease")
melt_tt$variable <- convert_names(melt_tt$variable, cbind(colnames(tt_df)[1:4], c("Conc.", "Hazard", "AUC", "OR")))

the_plot <- ggplot(melt_tt, aes(x = value, y = disease, )) + geom_point() + 
  facet_grid(cols = vars(variable), scales = "free") +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Statistic Difference")
plot(the_plot)
ggsave(paste0("reg_plots/reg_coef.tt.png"), the_plot, height = 6, width = 12)

tt_df <- data.frame(meta_results[[8]])
colnames(tt_df) <- c("diff_conc", "diff_survfit", "diff_auc", "diff_or")
tt_df$disease <- convert_names(meta_results[[1]])
tt_df <- tt_df[,c(5,1,2,3,4)]
tt_df[,2] <- signif(tt_df[,2], 3)
tt_df[,3] <- signif(tt_df[,3], 3)
tt_df[,4] <- signif(tt_df[,4], 3)
tt_df[,5] <- signif(tt_df[,5], 3)
write.table(tt_df, "supp_tables/train-test.txt", row.names = F, col.names = T, sep = "\t", quote = F)
