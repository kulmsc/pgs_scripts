
args <- commandArgs(trailingOnly=TRUE)

d <- args[1]
top_choice <- args[2]

#must make sure the number of requested top ranks SNPs is less than the total number of available SNPs
ss <- read.table(paste0(d, "/specific_ss"), stringsAsFactors=F)

if(top_choice < nrow(ss)){

  #construct a matrix with ncols = number of SNPs, nrows = 100 (sample size for each SNP)
  norm_samples <- apply(ss, 1, function(x) rnorm(100, mean = as.numeric(x[7]), sd = as.numeric(x[6])))

  #for each snp generate a rank of the SNP, note that the dimensions are now swapped from above
  rank_samples <- apply(norm_samples, 1, function(x) rank(x * -1))

  #calculate the number of time for each SNP the rank is less than the desired top number of SNPs
  prob_ranks <- apply(rank_samples, 1, function(x) sum(x < top_choice)/nrow(rank_samples))

  #multiple the normal beta values by the winner's curse probability to get an updated beta value
  ss[,7] <- ss[,7] * prob_ranks

  write.table(ss, paste0(d, "/adjust_ss"), col.names = F, row.names = F, sep = '\t', quote = F)
}
