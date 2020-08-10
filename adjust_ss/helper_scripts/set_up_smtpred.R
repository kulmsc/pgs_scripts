args <- commandArgs(trailingOnly=TRUE)

author <- args[1]
d <- args[2]
toprank <- as.numeric(args[3])

gencorr <- read.table("~/athena/doc_score/raw_ss/genetic_correlations", stringsAsFactors=F)
metastats <- read.table("~/athena/doc_score/raw_ss/meta_stats", stringsAsFactors=F, sep = ",", header = T)

gencorr <- gencorr[gencorr[,1] == author,]
gencorr <- gencorr[abs(gencorr[,8]) < gencorr[,7] & !is.na(gencorr[,7]) & gencorr[,2] != author,]
gencorr <- gencorr[order(gencorr[,7], decreasing=T),]

if(nrow(gencorr) > toprank){
  gencorr <- gencorr[1:toprank,]
}

metastats <- metastats[tolower(metastats[,1]) %in% c(author, gencorr[,2]),]
sampsize <- cbind(tolower(metastats$author), 4/(1/metastats$cases + 1/metastats$controls))
herit <- cbind(tolower(metastats$author), metastats$h2)
ss_names <- paste0(tolower(metastats$author), ".ss")

write.table(sampsize, paste0(d, "/samp_size"), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(herit, paste0(d, "/herit"), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(gencorr[,c(1,2,7)], paste0(d, "/rg"), row.names = F, col.names = F, quote = F, sep = '\t')
write.table(ss_names, paste0(d, "/ss_files"), row.names = F, col.names = F, quote = F, sep = '\t')
