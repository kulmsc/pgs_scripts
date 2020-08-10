#arguments are: d variable, i number, low author, chr
args <- commandArgs(trailingOnly=TRUE)
low_author <- tolower(args[3])

smtpred_ss <- read.table(paste0(args[1], "/multi_trait.beta"), stringsAsFactors=F, header = T)
raw_ss <- read.table(paste0(args[1], "/", args[3], ".ss"), stringsAsFactors=F, header=T)

raw_ss <- raw_ss[raw_ss$RSID %in% smtpred_ss$snpid,]
smtpred_ss <- smtpred_ss[order(smtpred_ss$snpid)[rank(raw_ss$RSID)],]

raw_ss$BETA <- smtpred_ss[,3]
write.table(raw_ss, paste0("~/athena/doc_score/mod_sets/", tools::toTitleCase(args[3]), "/", args[3], ".", args[4],  ".smtpred.", args[2] ,".ss"), row.names = F, col.names = T, sep = '\t', quote = F)
