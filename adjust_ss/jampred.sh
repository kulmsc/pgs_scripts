chr=$1
author=$2
dir=$3
d=comp_zone/dir${dir}
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`

#can limit to 5,000 samples and R2 < 0.95

head -5000 geno_files/${low_author}.${chr}.fam | cut -f1 > ${d}/subset_inds
plink --bfile geno_files/${low_author}.${chr} --keep-fam ${d}/subset_inds --indep-pairwise 500 100 0.95 --out ${d}/snps
plink --bfile geno_files/${low_author}.${chr} --keep-fam ${d}/subset_inds --extract ${d}/snps.prune.in --make-bed --out ${d}/for_jampred

Rscript helper_scripts/jampred.R $chr $author $d
