cat common_files/list_authors | while read author;do
#author=Tsoi
  low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
  ldsc --h2 ${author}/${low_author}.munged.sumstats.gz --ref-ld-chr ~/athena/refs/eur_w_ld_chr/ --w-ld-chr ~/athena/refs/eur_w_ld_chr/ --out ${author}/${low_author}.h2

  Rscript ~/Programs/HDL/HDL.run.R \
    gwas.df=${author}/${low_author}.munged.sumstats.gz \
    LD.path=~/athena/refs/UKB_imputed_SVD_eigen99_extraction \
    output.file=${author}/${low_author}.HDL.h2.Rout
done
