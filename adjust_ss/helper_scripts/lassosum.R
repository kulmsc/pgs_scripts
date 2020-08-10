#args are: d, chr, low_author
args <- commandArgs(trailingOnly=TRUE)

ss <- read.table(paste0("temp_files/ss.", args[3], ".", args[2]), stringsAsFactors=F, header=T)
ref.bfile <- paste0("geno_files/", args[3], ".", args[2])

cor <- p2cor(p = ss$P, n = ss$ESS[1], sign=ss$BETA)

out <- lassosum.pipeline(cor=cor, chr=ss$CHR, pos=ss$BP, 
                         A1=ss$A1, A2=ss$A2, # A2 is not required but advised
                         ref.bfile=ref.bfile, LDblocks = "EUR.hg19", sample = 5000,
                         s = c(0.1, 0.4, 0.8), lambda = c(0.002, 0.004, 0.007, 0.1))

new_ss <- ss[ss$BP %in% out$sumstats[,2],]

k <- 1
for(i in 1:3){
  for(j in 1:4){
    new_ss$BETA <- out$beta[[i]][,j]
    write.table(new_ss, paste0("~/athena/doc_score/mod_sets/", tools::toTitleCase(args[3]), "/", args[3], ".", args[2],  ".lassosum.", k ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
    k <- k + 1

  }
} 
