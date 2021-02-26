library(stringr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(viridis)
library(pROC)
theme_set(theme_cowplot())

author <- "Bentham"
all_authors <- unique(str_split(list.files("nonpred_results/per_score/"), fixed("."), simplify = T)[,1])
all_authors <- all_authors[-which(all_authors == "xie")]

for(author in all_authors){
  print(author)

nonpred_res <- readRDS(paste0("nonpred_results/per_score/", tolower(author), ".res.many.RDS"))
score_names <- nonpred_res[["extra"]][[1]]
method_names <- str_split(score_names, fixed("."), simplify = T)[,3]
best_in_method <- read.table(paste0("tune_results/", tolower(author), ".methods.ss"))

convert_names <- function(x, the_dict = method_dict){
  y <- rep(NA, length(x))
  for(i in 1:length(x)){
    if(x[i] %in% the_dict[,1]){
      y[i] <- the_dict[the_dict[,1] == x[i], 2]
    }
  }
  return(y)
}

method_dict <- read.table("local_info/method_names", stringsAsFactors = F)
disease_dict <- read.table("local_info/disease_names", stringsAsFactors = F, sep = "\t")

neat_method_names <- convert_names(method_names)
neat_score_names <- paste0(neat_method_names, "-", str_split(score_names, fixed("."), simplify = T)[,2])

