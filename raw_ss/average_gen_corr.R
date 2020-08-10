
gencorr <- read.table("genetic_correlations", stringsAsFactors=F, fill = T)
gencorr$avg <- NA
gencorr$se <- NA

for(i in 1:nrow(gencorr)){

  if(is.finite(gencorr[i,3]) & is.finite(gencorr[i,5])){
    weights <- c(1/gencorr[i,4], 1/gencorr[i,6]) / sum(c(1/gencorr[i,4], 1/gencorr[i,6]))
    gencorr$avg[i] <- sum(c(gencorr[i,3], gencorr[i,5]) * weights)
    gencorr$se[i] <- sum(c(gencorr[i,4], gencorr[i,6]) * weights)

  } else if(is.finite(gencorr[i,3])){
    gencorr$avg[i] <- gencorr[i,3]
    gencorr$se[i] <- gencorr[i,4]

  } else if(is.finite(gencorr[i,5])){
    gencorr$avg[i] <- gencorr[i,5]
    gencorr$se[i] <- gencorr[i,6]
  } 

}

write.table(gencorr, "genetic_correlations", row.names = F, col.names = F, sep = '\t', quote = F)
