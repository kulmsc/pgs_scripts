library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(viridis)
theme_set(theme_cowplot())


method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
anno_dict <- read.table("local_info/anno_names", stringsAsFactors = F, sep = "\t")
author_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = "\t")

convert_names <- function(x, the_dict = method_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}



########################################################################################################
#                     FUNC ANNO WEIGHT  ###############################################################
########################################################################################################


trait_name <- convert_names(read.table("anno_scores/anno_names")[,1], anno_dict)
all_author <- unique(str_split(list.files("anno_scores/"), fixed("."), simplify = T)[,1])[-1]
all_author <- all_author[all_author != "Xie"]

total_anno_score <- matrix(0, nrow = 15, ncol = 44)
total_score_name <- c("Clump", "DBLSMM", "Double-Weight", "JAMPred", "lassosum", "LDpred", "LDpred2",
                      "prsCS", "SBayesR", "SBLUP", "SMTpred", "Tweedie", "WC-2D", "WC-Lasso", "WC-Likelihood")
goes_into <- rep(0, length(total_score_name))
best_anno_score <- list()
ii <- 1

for(author in all_author){
    print(author)
    score_name <- read.table(paste0("tune_results/", tolower(author), ".methods.ss"))
    score_name <- convert_names(str_split(score_name[,1], fixed("."), simplify = T)[,3])
    best_score <- read.table(paste0("tune_results/", tolower(author), ".best.ss"), stringsAsFactors = F)
    best_score <- convert_names(strsplit(best_score[3,2], ".", fixed=T)[[1]][3])
    print(length(score_name))
    
    anno_score <- read.table(paste0("anno_scores/", author, ".anno.score"), stringsAsFactors = F)
    anno_total <- read.table(paste0("anno_scores/", author, ".anno.total"), stringsAsFactors = F)
    anno_score <- anno_score/anno_total
    
    for(i in 1:nrow(anno_score)){
      anno_score[i,] <- (anno_score[i,] - min(anno_score[i,], na.rm = T))/(max(anno_score[i,], na.rm = T) - min(anno_score[i,], na.rm = T))
      anno_score[i,] <- anno_score[i,]/sum(anno_score[i,], na.rm = T)
    }
    
    for(i in 1:nrow(anno_score)){
      goes_into[which(total_score_name == score_name[i])] <- goes_into[which(total_score_name == score_name[i])] + 1
      for(j in 1:ncol(anno_score)){
        if(!is.na(anno_score[i,j])){
          total_anno_score[which(total_score_name == score_name[i]), j] <- total_anno_score[which(total_score_name == score_name[i]), j] +
            anno_score[i,j]
        }
      }
    }
    
    print(dim(anno_score))
    colnames(anno_score) <- trait_name
    anno_score$method <- score_name

    anno_score <- melt(anno_score, id.vars = "method")
    
    the_plot <- ggplot(anno_score, aes(method, variable, fill = value)) + geom_raster() +
      scale_fill_viridis(begin = 0, end = 1) + 
      labs(x = "", y = "", fill = "Func.\nAnno.\nWeight") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plot(the_plot)
    ggsave(paste0("anno_plots/anno.", author, ".png"), the_plot, height = 12, width = 10)

    anno_score$author <- author
    best_anno_score[[ii]] <- anno_score[anno_score$method == best_score,]
    ii <- ii + 1
}

########################################################################################################
#              SUMMARY PLOTS
########################################################################################################

fa_order <- as.character(best_anno_score[[1]]$variable)
arrange_best_anno_score <- do.call("cbind", lapply(best_anno_score, function(x) x[,3,drop=F]))
x <- apply(arrange_best_anno_score, 1, function(x) mean(x, na.rm = T))

for_later <- as.character(best_anno_score[[1]]$variable)
or1 <-  hclust(dist(arrange_best_anno_score))$order
or2 <- hclust(dist(t(arrange_best_anno_score)))$order

for_supp <- data.frame(convert_names(all_author, author_dict), t(arrange_best_anno_score))
colnames(for_supp) <- c("Disease", trait_name)
for(i in 2:ncol(for_supp)){
  for_supp[,i] <- signif(for_supp[,i], 3)
}
write.table(for_supp[,1:8], "supp_tables/anno.1.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,9:15)], "supp_tables/anno.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,16:22)], "supp_tables/anno.3.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,23:29)], "supp_tables/anno.4.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,30:36)], "supp_tables/anno.5.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,37:42)], "supp_tables/anno.6.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,43:45)], "supp_tables/anno.7.txt", row.names = F, col.names = T, quote = F, sep = "\t")


best_anno_score <- do.call("rbind", best_anno_score)
best_anno_score$disease <- convert_names(best_anno_score$author, author_dict)
best_anno_score$disease <- factor(best_anno_score$disease, levels = unique(best_anno_score$disease)[or2])
best_anno_score$variable <- factor(best_anno_score$variable, levels = for_later[or1])

the_plot <- ggplot(best_anno_score, aes(disease, variable, fill = value)) + geom_raster() +
  scale_fill_viridis(trans = scales::pseudo_log_trans(sigma = 0.02)) +
  labs(x = "", y = "", fill = "Func.\nAnno.\nWeight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(the_plot)
ggsave(paste0("anno_plots/anno.disease.png"), the_plot, height = 12, width = 10)

exclude_vars <- c("Dnase I Hypersensitive-Fetal", "TFBS", "Other Cell Type",
                  "Weak Enhancer", "Promoter Flanking", "GI Cell Type", "H3K27ac",
                  "Connective Bone Cell Type", "Skeletal/Muscle Cell Type",
                  "Adrenal/Pancreas Cell Type", "Dispersed Gene Family")
plot_df <- best_anno_score[!(best_anno_score$variable %in% exclude_vars),]
the_plot <- ggplot(plot_df, aes(disease, variable, fill = value)) + geom_raster() +
  scale_fill_viridis(trans = scales::pseudo_log_trans(sigma = 0.02)) +
  labs(x = "", y = "", fill = "Func. Anno. Weight:") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "top",
        legend.key.width = unit(1.9, 'cm'),
        legend.spacing.x = unit(0.6, 'cm')) 
  
plot(the_plot)
ggsave(paste0("anno_plots/anno.SMALL.disease.png"), the_plot, height = 10, width = 8.5)

#Methods
for(i in 1:nrow(total_anno_score)){
  total_anno_score[i,] <- total_anno_score[i,]/goes_into[i]
}
total_anno_score <- data.frame(total_anno_score)
colnames(total_anno_score) <- trait_name
total_anno_score$method <- total_score_name
or1 <-  hclust(dist(total_anno_score[,-ncol(total_anno_score)]))$order
or2 <- hclust(dist(t(total_anno_score[,-ncol(total_anno_score)])))$order

anno_score <- melt(total_anno_score, id.vars = "method")

anno_score$method <- factor(anno_score$method, levels = total_score_name[or1])
anno_score$variable <- factor(anno_score$variable, levels = as.character(unique(anno_score$variable))[or2])


the_plot <- ggplot(anno_score, aes(method, variable, fill = value)) + geom_raster() +
  scale_fill_viridis() + labs(x = "", y = "", fill = "Func.\nAnno.\nWeight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(the_plot)
ggsave(paste0("anno_plots/anno.total_meth.png"), the_plot, height = 12, width = 10)



########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################



########################################################################################################
#                     DELETERIOUS SCORE ###############################################################
########################################################################################################

trait_name <- convert_names(read.table("anno_scores/anno_names")[,1], anno_dict)
all_author <- unique(str_split(list.files("anno_scores/"), fixed("."), simplify = T)[,1])[-1]
all_author <- all_author[all_author != "Xie"]

total_anno_score <- matrix(0, nrow = 15, ncol = 10)
total_score_name <- c("Clump", "DBLSMM", "Double-Weight", "JAMPred", "lassosum", "LDpred", "LDpred2",
                      "prsCS", "SBayesR", "SBLUP", "SMTpred", "Tweedie", "WC-2D", "WC-Lasso", "WC-Likelihood")
score_type <- c('SIFT', 'Polyphen2-HDIV', 'Polyphen2-HVAR', 'LRT', 'MutationTaster',
                'FATHMM', 'PROVEAN', 'VEST4', 'MVP', 'ClinPred')
goes_into <- rep(0, length(total_score_name))
best_anno_score <- list()
ii <- 1

for(author in all_author){
    
    anno_score <- read.table(paste0("anno_scores/", author, ".dbnsfp.score"))
    anno_total <- read.table(paste0("anno_scores/", author, ".dbnsfp.len"))
    anno_score <- anno_score/anno_total[,rep(1,10)]
    anno_score[is.na(anno_score)] <- 0
    
    for(i in 1:nrow(anno_score)){
      if(sum(anno_score[i,]) > 0){
        anno_score[i,] <- (anno_score[i,] - min(anno_score[i,]))/(max(anno_score[i,]) - min(anno_score[i,]))
        anno_score[i,] <- anno_score[i,]/sum(anno_score[i,])
      }
    }
    
    
    for(i in 1:nrow(anno_score)){
      goes_into[which(total_score_name == total_score_name[i])] <- goes_into[which(total_score_name == total_score_name[i])] + 1
      for(j in 1:ncol(anno_score)){
        if(!is.na(anno_score[i,j])){
          total_anno_score[which(total_score_name == total_score_name[i]), j] <- 
            total_anno_score[which(total_score_name == total_score_name[i]), j] +
            anno_score[i,j]
        }
      }
    }
    
    
    score_name <- read.table(paste0("tune_results/", tolower(author), ".methods.ss"))
    score_name <- convert_names(str_split(score_name[,1], fixed("."), simplify = T)[,3])
    
    colnames(anno_score) <- score_type
    anno_score$method <- score_name
    
    anno_score <- melt(anno_score, id.vars = "method")
    
    the_plot <- ggplot(anno_score, aes(method, variable, fill = value)) + geom_raster() +
      scale_fill_viridis() + labs(x = "", y = "", fill = "Func.\nAnno.\nWeight") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plot(the_plot)
    ggsave(paste0("anno_plots/delet.", author, ".png"), the_plot, height = 12, width = 10)

    anno_score$author <- author
    best_score <- read.table(paste0("tune_results/", tolower(author), ".best.ss"), stringsAsFactors = F)
    best_score <- convert_names(strsplit(best_score[3,2], ".", fixed=T)[[1]][3])
    best_anno_score[[ii]] <- anno_score[anno_score$method == best_score,]
    ii <- ii + 1
}


########################################################################################################
#              SUMMARY PLOTS
########################################################################################################

arrange_best_anno_score <- do.call("cbind", lapply(best_anno_score, function(x) x[,3,drop=F]))
for_later <- as.character(best_anno_score[[1]]$variable)
or1 <-  hclust(dist(arrange_best_anno_score))$order
or2 <- hclust(dist(t(arrange_best_anno_score)))$order
best_anno_score <- do.call("rbind", best_anno_score)
best_anno_score$disease <- convert_names(best_anno_score$author, author_dict)
best_anno_score$disease <- factor(best_anno_score$disease, levels = unique(best_anno_score$disease)[or2])
best_anno_score$variable <- factor(best_anno_score$variable, levels = for_later[or1])

for_supp <- data.frame(convert_names(all_author, author_dict), t(arrange_best_anno_score))
colnames(for_supp) <- c("Disease", score_type)
for(i in 2:ncol(for_supp)){
  for_supp[,i] <- signif(for_supp[,i], 3)
}
write.table(for_supp[,1:5], "supp_tables/delet.1.txt", row.names = F, col.names = T, quote = F, sep = "\t")
write.table(for_supp[,c(1,6:10)], "supp_tables/delet.2.txt", row.names = F, col.names = T, quote = F, sep = "\t")

the_plot <- ggplot(best_anno_score, aes(disease, variable, fill = value)) + geom_raster() +
  scale_fill_viridis() + labs(x = "", y = "", fill = "Func.\nAnno.\nWeight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(the_plot)
ggsave(paste0("anno_plots/delet.disease.png"), the_plot, height = 12, width = 10)

#Methods
for(i in 1:nrow(total_anno_score)){
  total_anno_score[i,] <- total_anno_score[i,]/goes_into[i]
}
total_anno_score <- data.frame(total_anno_score)
colnames(total_anno_score) <- score_type
total_anno_score$method <- total_score_name
or1 <-  hclust(dist(total_anno_score[,-ncol(total_anno_score)]))$order
or2 <- hclust(dist(t(total_anno_score[,-ncol(total_anno_score)])))$order

anno_score <- melt(total_anno_score, id.vars = "method")

anno_score$method <- factor(anno_score$method, levels = total_score_name[or1])
anno_score$variable <- factor(anno_score$variable, levels = score_type[or2])

the_plot <- ggplot(anno_score, aes(method, variable, fill = value)) + geom_raster() +
  scale_fill_viridis() + labs(x = "", y = "", fill = "Func.\nAnno.\nWeight") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(the_plot)
ggsave(paste0("anno_plots/delet.total_meth.png"), the_plot, height = 12, width = 10)
