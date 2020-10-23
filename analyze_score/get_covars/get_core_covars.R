
covars <- read.table("into_base_covars", stringsAsFactors=F, header=T, sep = ",")

dob <- as.Date(paste0(covars[,3], "/", covars[,4], "/", "15"))
age <- as.numeric((Sys.Date() - dob)/365)

covars <- data.frame(covars[,1:2], age, covars[,5:ncol(covars)])
covars <- covars[,-ncol(covars)]
colnames(covars) <- c("eid", "sex", "age", paste0("PC", 1:10))

saveRDS(covars, "base_covars.RDS")
write.table(covars, "base_covars", col.names = T, row.names = F, sep = '\t', quote = F)
