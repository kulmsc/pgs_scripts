
args <- commandArgs(trailingOnly=TRUE)

chr <- args[1]
author <- args[2]
d <- args[3]
i <- args[4]
lowauthor <- tolower(author)

sblup <- read.table(paste0(d,"/out.sblup.cojo"), stringsAsFactors=F)
ss <- read.table(paste0("temp_files/ss.", lowauthor, ".", chr), stringsAsFactors=F, header=T)

ss <- ss[ss$RSID %in% sblup[,1],]
sblup <- sblup[order(sblup[,1])[rank(ss$RSID)],]
ss$BETA <- sblup[,4]

write.table(ss, paste0("../mod_sets/", author, "/", lowauthor, ".", chr, ".sblup.", i, ".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
