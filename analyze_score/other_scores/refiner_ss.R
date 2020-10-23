library(vroom)
library(stringr)
library(bigsnpr)

options(warn=2)

#Read in the UKBB summary statistics
impute <- read.table("../../raw_ss/common_files/impute_rsids", header = F, stringsAsFactors = F)
colnames(impute) <- c("LOC", "SNP", "POS", "A1", "A2", "MAF", "AX", "INFO")

#Make sure the UKBB summary statistics have alleles that match ACGT and are not ambigous
impute <- impute[nchar(impute$A1) == 1 & nchar(impute$A2) == 1,]
impute <- impute[impute$A1 %in% c("A", "C", "G", "T") & impute$A2 %in% c("A", "C", "G", "T"),]
impute <- impute[!((impute$A1 == "A" & impute$A2 == "T") |
             (impute$A1 == "T" & impute$A2 == "A") |
             (impute$A1 == "G" & impute$A2 == "C") |
             (impute$A1 == "C" & impute$A2 == "G")),]

#Read in other UKBB data that comes with chromosomes so I can give a chromosome to each SNP
impute$CHR <- "X"
for(i in 1:22){
  temp <- readRDS(paste0("~/athena/ukbiobank/qc/imputed/chr", i, ".RDS"))
  impute$CHR[impute$SNP %in% temp[,2]] <- i
}
impute <- impute[impute$CHR %in% 1:22,]
impute$SUPERPOS <- paste0(impute$CHR, "_", impute$POS)

#Now iterating through each of the score summary statistics
###
for(filename in list.files("raw_other_ss/")){
  #read the summary statistic in
  ss <- read.table(paste0("raw_other_ss/", filename), stringsAsFactors=F, header=T, sep='\t')

  #require both a column for major and minor allele, if only one exists refining is not possible
  if("effect_allele" %in% colnames(ss) & "reference_allele" %in% colnames(ss)){

    #limit SNPs to those with ACGT
    ss$effect_allele <- toupper(ss$effect_allele)
    ss$reference_allele <- toupper(ss$reference_allele)
    ss <- ss[nchar(ss$effect_allele) == 1 & nchar(ss$reference_allele) == 1,]
    ss <- ss[ss$effect_allele %in% c("A", "C", "G", "T") & ss$reference_allele %in% c("A", "C", "G", "T"),]

    #remove ambigous SNPs
    ss <- ss[!((ss$effect_allele == "A" & ss$reference_allele == "T") |
               (ss$effect_allele == "T" & ss$reference_allele == "A") |
               (ss$effect_allele == "G" & ss$reference_allele == "C") |
               (ss$effect_allele == "C" & ss$reference_allele == "G")),]

    #if rsID in the ss use it to match and sort to the UKBB
    #if not then use the chromosome and position (I already checked that all scores' ref genomes match UKBB)
    if("rsID" %in% colnames(ss)){
      ss <- ss[!is.na(ss$rsID),]
      ss$rsID <- tolower(ss$rsID)
      ss <- ss[ss$rsID %in% impute$SNP,]
      sub_impute <- impute[impute$SNP %in% ss$rsID,]
      ss <- ss[order(ss$rsID)[rank(sub_impute$SNP)],]
    } else {
      ss$SUPERPOS <- paste0(ss$chr_name, "_", ss$chr_position)
      ss <- ss[ss$SUPERPOS %in% impute$SUPERPOS,]
      sub_impute <- impute[impute$SUPERPOS %in% ss$SUPERPOS,]
      ss <- ss[order(ss$SUPERPOS)[rank(sub_impute$SUPERPOS)],]
      ss$rsID <- sub_impute$SNP    
    }

    #Change the column names and then use the snp_match function from bigsnpr to flip alleles
    colnames(sub_impute) <- c("loc", "rsid", "pos", "a0", "a1", "maf", "ax", "info", "chr", "superpos")
    colnames(ss)[colnames(ss) == "effect_allele"] <- "a0"
    colnames(ss)[colnames(ss) == "reference_allele"] <- "a1"
    colnames(ss)[colnames(ss) == "rsID"] <- "rsid"
    colnames(ss)[colnames(ss) == "effect_weight"] <- "beta"
    ss$pos <- sub_impute$pos
    ss$chr <- sub_impute$chr
    ss <- snp_match(ss, sub_impute)

    #write out the refined summary statistic
    out_ss <- data.frame(chr = ss$chr, rsid = ss$rsid.ss, A1 = ss$a0, beta = ss$beta)
    new_name <- strsplit(filename, ".", fixed = T)[[1]][1]
    write.table(out_ss, paste0("refine_other_ss/", new_name, ".new.ss"), col.names = T, row.names = F, sep = '\t', quote=F)
  }
}
