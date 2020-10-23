
ss_name=$1
chr=$2
i=$3

echo $ss_name $author $method $chr $i

#Just logistics around controlling parallelism
go_num=`head -1 temp_files/poss_go`
grep -v -w $go_num temp_files/poss_go > temp_files/temp_poss
mv temp_files/temp_poss temp_files/poss_go

if [ ! -e other_small_score_files/score.${ss_name}.${chr}.profile.zst ];then

  len=`cat refine_other_ss/${ss_name} | awk -v var="$chr" '$1 == var {print $0}' | wc -l`
  if [ $len != 0 ];then
    cat refine_other_ss/${ss_name} | awk -v var="$chr" '$1 == var {print $0}' > temp_files/ss.${i}

    cat temp_files/ss.${i} | cut -f2 > temp_files/rsids.${i}

    bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chr}.bgen -incl-rsids temp_files/rsids.${i} > temp_files/temp.${i}.bgen

    plink2_new --memory 12000 --threads 12 --bgen temp_files/temp.${i}.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chr}.sample --keep-fam temp_files/nonbrit_eid --make-bed --out temp_files/geno.${i}

    plink --memory 12000 --threads 12 --bfile temp_files/geno.${i} --keep-allele-order --score temp_files/ss.${i} 2 3 4 sum --out other_small_score_files/score.${ss_name}.${chr}

    zstd --rm other_small_score_files/score.${ss_name}.${chr}.profile
  fi
fi

  rm temp_files/ss.${i}
  rm temp_files/rsids.${i}
  rm temp_files/temp.${i}.bgen
  rm temp_files/geno.${i}.*


#Just logistics around controlling parallelism
echo $go_num >> temp_files/poss_go

