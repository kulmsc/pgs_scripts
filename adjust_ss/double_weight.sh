chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

plink --memory 4000 --threads 1 --bfile geno_files/${low_author}.${chr} --clump temp_files/ss.${low_author}.${chr} --clump-snp-field RSID --clump-p1 0.5 --clump-r2 0.05 --out ${d}/out
sed -e 's/ [ ]*/\t/g' ${d}/out.clumped | sed '/^\s*$/d' | cut -f4 | tail -n +2 > ${d}/done_rsids
cat temp_files/ss.${low_author}.${chr} | fgrep -w -f ${d}/done_rsids > ${d}/specific_ss

i=1

cat all_specs/double_param_specs | while read spec;do
  if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.double.${i}.ss ]; then
    Rscript helper_scripts/double_weight.R $d $spec
    if [ -e ${d}/adjust_ss ]; then
      mv ${d}/adjust_ss ~/athena/doc_score/mod_sets/${author}/${low_author}.${chr}.double.${i}.ss
    fi
  fi
  let i=i+1
done
