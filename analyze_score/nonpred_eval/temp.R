list_all_best_names <- list()
list_all_best_methods <- list()

all_author <- unlist(lapply(strsplit(list.files("../tune_score/tune_results/", "res"), "_"), function(x) x[1]))
all_author <- all_author[all_author != "Xie"]
i <- 1
for(author in all_author){

  list_all_best_names[[i]] <- read.table(paste0("../tune_score/tune_results/", tolower(author), ".best.ss"), stringsAsFactors=F)
  list_all_best_methods[[i]] <- read.table(paste0("../tune_score/tune_results/", tolower(author), ".methods.ss"), stringsAsFactors=F)
  i <- i + 1
}

abs_best_names <- unlist(lapply(list_all_best_names, function(x) x[3,2]))
saveRDS(abs_best_names, "best_names.RDS")
