library(stringr)

all_files <- list.files("from_athena/", "RDS")

out_scores <- replicate(3,list(NULL))
all_authors <- c()

for(f in all_files){
 all_scores <- readRDS(paste0("from_athena/", f))
 current_authors <- str_split(colnames(all_scores), fixed("."), simplify = T)[,1]
 all_authors <- unique(c(all_authors, current_authors))
 for(author in unique(current_authors)){
  out_scores[[which(all_authors == author)]] <- cbind(out_scores[[which(all_authors == author)]], 
                                                      all_scores[,current_authors == author])
 }
}


for(author in all_authors){
 saveRDS(out_scores[[which(all_authors == author)]], paste0("from_athena/", author, ".tune.RDS"))
}

