chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

num_snps=`cat ../raw_ss/meta_stats | fgrep $author | cut -f7 -d','`
h2=`cat ../raw_ss/meta_stats | fgrep $author | cut -f12 -d','`
samp_size=`cat temp_files/ss.${low_author}.${chr} | head -2 | tail -n +2 | cut -f9`

i=1
cat all_specs/dbslmm_param_specs | while read spec;do

  if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.dbslmm.${i}.ss ]; then

    plink --bfile geno_files/${low_author}.${chr} --freq --out temp_files/${low_author}.${chr}
    Rscript helper_scripts/add_gemma.R $author $chr

    dbslmm_bash=~/Programs/DBSLMM/dbslmm
    dbslmm_r=~/Programs/DBSLMM/software/DBSLMM.R
    ss=temp_files/ss.${low_author}.${chr}.gemma
    ref=geno_files/${low_author}.${chr}
    blockf=~/athena/refs/ldsplits/chr${chr}.bed 
    
    Rscript ${dbslmm_r} --summary $ss --outPath ${d}/ --plink ~/bin/plink --dbslmm ${dbslmm_bash} --ref ${ref} --n ${samp_size} --nsnp ${num_snps} --block ${blockf} --h2 ${h2}

  fi

  let i=i+1
done
