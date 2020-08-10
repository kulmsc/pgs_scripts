library(bigsnpr)
library(R2BGLiMS)

args <- commandArgs(trailingOnly=TRUE)

chr <- args[1]
author <- args[2]
d <- args[3]

if(!file.exists("comp_zone/dir13/for_jampred.rds")){
  snp_readBed("comp_zone/dir13/for_jampred.bed")
}
obj.bigSNP <- snp_attach(paste0("comp_zone/dir13/for_jampred.rds"))
ped <- big_copy(snp_fastImputeSimple(obj.bigSNP$genotypes, "mean0"), type = "integer")
ped <- ped$bm()
options(bigmemory.allow.dimnames=TRUE)
colnames(ped) <- obj.bigSNP$map$marker.ID
rm(obj.bigSNP)

ss <- read.table("temp_files/ss.bentham.22", stringsAsFactors=F, header=T)
ss <- ss[ss$RSID %in% colnames(ped),]
marg_beta <- ss$BETA
names(marg_beta) <- ss$RSID
marg_se <- ss$SE
names(marg_se) <- ss$RSID

meta_stats <- read.table("~/athena/doc_score/raw_ss/meta_stats", sep = ",", stringsAsFactors=F, header = T)
meta_line <- meta_stats[meta_stats[,1] == tools::toTitleCase(author),]


jampred.res.bin <- JAMPred(
 marginal.betas = marg_beta,
 n.training = as.numeric(meta_line$sampe_size),
 marginal.logor.ses = marg_se, # Only necessary for a binary trait
 p.cases.training = meta_line$cases/meta_line$sampe_size, # Only necessary for a binary trait
 ref.geno = ped,
 total.snps.genome.wide = meta_line$snps, # Total SNPs across all chromosomes
 n.mil = 0.2,
 seed = 1 # For re-producibility. If not set a random seed is used
)
