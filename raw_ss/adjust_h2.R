info <- read.csv("meta_stats", stringsAsFactors=F, header=T)

info$hdl_h2[grepl("was", info$hdl_h2)] <- 0
info$hdl_h2se[grepl("estimated", info$hdl_h2se)] <- 0

info$hdl_h2 <- as.numeric(info$hdl_h2)
info$hdl_h2se <- as.numeric(info$hdl_h2se)

weights <-  cbind(1/info$ldsc_h2_se, 1/info$hdl_h2se)/rowSums(cbind(1/info$ldsc_h2_se, 1/info$hdl_h2se))
new_h2 <- rowSums(cbind(info$ldsc_h2, info$hdl_h2) * weights)

final_h2 <- rep(0, nrow(info))
final_h2[info$hdl_h2 == 0] <- info$ldsc_h2[info$hdl_h2 == 0]
final_h2[info$ldsc_h2 < 0 | info$ldsc_h2 > 1] <- info$hdl_h2[info$ldsc_h2 < 0 | info$ldsc_h2 > 1]
final_h2[final_h2 == 0] <- new_h2[final_h2 == 0]
final_h2[is.nan(final_h2)] <- 0.01

info$h2 <- final_h2
write.table(info, "meta_stats", row.names = F, col.names = T, quote = F, sep = ',')
