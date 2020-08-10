#double for loop to get pairwise correlations
#always remember to load the right python environment

cat common_files/temp_authors | while read author;do
  cat common_files/list_authors | while read sec_author;do

  low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
  o_author=`echo "$sec_author" | tr '[:upper:]' '[:lower:]'`

  ldsc --rg ${author}/${low_author}.munged.sumstats.gz,${sec_author}/${o_author}.munged.sumstats.gz --ref-ld-chr ~/athena/refs/eur_w_ld_chr/ --w-ld-chr ~/athena/refs/eur_w_ld_chr/ --out gen_corr/ldsc/${low_author}.${o_author}.corr

  done
done
