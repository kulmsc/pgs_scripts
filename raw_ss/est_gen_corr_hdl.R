library(HDL)

all_authors <- read.table("common_files/list_authors", stringsAsFactors = F)
gwas_files <- list()
LD.path <- "~/athena/refs/UKB_imputed_SVD_eigen99_extraction/"

# read in all of the authors once up top
i <- 1
for(author in all_authors[,1]){
  low_author <- tolower(author)
  gwas_files[[i]] <- read.table(paste0(author, "/", low_author, ".munged.sumstats.gz"), stringsAsFactors=F, header=T)
  i <- i + 1
}

# create a doule for loop to get pairwise correlations
for(i in 12:length(gwas_files)){
  for(j in 1:length(gwas_files)){
    if(i != j & (!(i == 12 & j %in% 1:24))){ #I was having some trouble completing some calculations so I would manually set this change
      res.HDL <- HDL.rg(gwas_files[[i]], gwas_files[[j]], LD.path)
      write.table(res.HDL$estimates.df, paste0("gen_corr/hdl/", tolower(all_authors[i,1]), ".", tolower(all_authors[j,1]), ".corr.log"), row.names = T, col.names = T, quote = F)
    }
  }
}
