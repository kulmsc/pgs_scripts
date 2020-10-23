library(data.table)

all_files <- list.files("small_score_files/", pattern = "profile")
split_files <- strsplit(all_files, ".", fixed = T)
all_author <- unlist(lapply(split_files, function(x) x[2]))
all_chr <- unlist(lapply(split_files, function(x) x[5]))

brit_eid <- read.table("temp_files/brit_eid", stringsAsFactors=F)
new_name <- unique(all_author)
all_score <- matrix(0, nrow = nrow(brit_eid), ncol = length(new_name))
colnames(all_score) <- new_name


i <- 1
j <- 1
for(i in 1:length(all_files)){

          system(paste0("zstd -d small_score_files/", all_files[i]))
          sub_score <- as.data.frame(fread(paste0("small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4))))
          system(paste0("rm small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4)))
          all_score[,new_name == all_author[i]] <- all_score[,new_name == all_author[i]] + sub_score$SCORESUM 
    
}

all_score <- all_score[,colSums(all_score) != 0]

next_file <- length(grep("RDS", list.files("final_scores", "all_score"))) + 1

write.table(all_score, paste0("final_scores/all_score.", next_file), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(all_score, paste0("final_scores/all_score.", next_file, ".RDS"))
write.table(colnames(all_score), paste0("final_scores/done_names.", next_file), row.names = F, col.names = F, quote = F)
system(paste0("gzip final_scores/all_score.", next_file))
