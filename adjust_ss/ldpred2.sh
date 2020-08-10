chr=$1
author=$2
dir=$3
d=comp_zone/dir${dir}
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`

head -10000 geno_files/${low_author}.${chr}.fam | cut -f1 > ${d}/subset_inds
plink --bfile geno_files/${low_author}.${chr} --keep-fam ${d}/subset_inds --extract ~/athena/refs/hapmap_from_ldpred2 --make-bed --out ${d}/for_ldpred2

Rscript helper_scripts/ldpred2.R $chr $author $d
