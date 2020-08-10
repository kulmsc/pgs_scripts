library(bigsnpr)
#source("helper_scripts/snp_as_genetic.R")

args <- commandArgs(trailingOnly=TRUE)

chr <- args[1]
author <- args[2]
d <- args[3]

lowauthor <- tolower(author)

#read in the genotypic data
if(!file.exists(paste0(d, "/for_ldpred2.rds"))){
  snp_readBed(paste0(d, "/for_ldpred2.bed"))
}
obj.bigSNP <- snp_attach(paste0(d, "/for_ldpred2.rds"))

#split up genotypic data as described
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection - 1
NCORES <- nb_cores()

#read in the summary stats
sumstats <- bigreadr::fread2(paste0("temp_files/ss.", lowauthor, ".", chr))
colnames(sumstats) <- c("chr", "pos", "rsid", "a0", "a1", "beta_se", "beta", "p", "n_eff")
sumatsts <- sumstats[sumstats$beta_se > 0,]

#run snp_match (likely can skip since I do this earlier)
map <- obj.bigSNP$map[-(2:3)]
names(map) <- c("chr", "pos", "a0", "a1")
info_snp <- snp_match(sumstats, map, join_by_pos = T)

#Calculate the LD correlation matrix
print("calc the corr matrix")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "temp_files", ncores = 3)
df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
corr0 <- snp_cor(G, ncores = NCORES, infos.pos = POS2, size = 3 / 1000)
corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))

#Get the heritability
print("getting h2")
h2 <- read.table("../raw_ss/meta_stats", stringsAsFactors=F, sep = ",", header=T)
h2 <- h2$h2[h2$author == author]

#Actually run ldpred2
print("run ldpred2 inf")
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2)

#set up parameters for the non-inf version
h2_seq <- round(h2 * c(0.7, 1, 1.4), 3)
p_seq <- signif(seq_log(1e-4, 1, length.out = 4), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

#run ldpred2 over all the parameters
print("run ldpred2 over grid")
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)

#write the results to mod_sets
print("writing the results")
ss <- info_snp[,1:9]
ss$beta <- beta_inf
write.table(ss, paste0("../mod_sets/", author, "/", lowauthor, ".", chr,  ".ldpred2.1.ss"), row.names = F, col.names = T, quote = F, sep = '\t')

for(i in 1:ncol(beta_grid)){
  ss$beta <- beta_grid[,i]
  write.table(ss, paste0("../mod_sets/", author, "/", lowauthor, ".", chr,  ".ldpred2.", i+1, ".ss"), row.names = F, col.names = T, quote = F, sep = '\t')
}
