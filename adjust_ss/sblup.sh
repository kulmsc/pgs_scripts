chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

num_snps=`cat ../raw_ss/meta_stats | fgrep $author | cut -f7 -d','`
h2=`cat ../raw_ss/meta_stats | fgrep $author | cut -f12 -d','`

i=1
cat all_specs/sblup_param_specs | while read spec;do
  coef=`echo $spec | cut -f1 -d' '`
  fam_num=`echo $spec | cut -f2 -d' '`
  use_h2=`echo "$h2 * $coef" | bc -l`
  lambda=`echo "$num_snps * (1/$use_h2 - 1)" | bc`

  if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.sblup.${i}.ss ]; then

    plink --bfile geno_files/${low_author}.${chr} --freq --out temp_files/${low_author}.${chr}
    Rscript helper_scripts/add_ma.R $author $chr

    head -$fam_num geno_files/${low_author}.${chr}.fam | cut -f1 > ${d}/subset_inds
    plink --bfile geno_files/${low_author}.${chr} --keep-fam ${d}/subset_inds --make-bed --out ${d}/for_sblup

    gcta64 --bfile ${d}/for_sblup --cojo-file temp_files/ss.${low_author}.${chr}.ma --cojo-sblup $lambda --cojo-wind 100 --thread-num 20 --out ${d}/out
    Rscript helper_scripts/reformat_sblup.R $chr $author $d $i
    rm ${d}/for_sblup*

  fi
  let i=i+1
done
