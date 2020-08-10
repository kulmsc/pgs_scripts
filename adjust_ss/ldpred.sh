chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}

prev_reftype="none"

i=1
cat all_specs/ldpred_param_specs | tail -n +2 | while read spec;do
  fval=`echo $spec | cut -f1 -d' '`
  reftype=`echo $spec | cut -f2 -d' '`

  echo fval $fval
  echo reftype $reftype

  num_snps=`cat ../raw_ss/meta_stats | fgrep $author | cut -f7 -d','`
  h2=`cat ../raw_ss/meta_stats | fgrep $author | cut -f12 -d','`
  rad=`echo $num_snps/3000 | bc`

  if [ $reftype != $prev_reftype ];then
    echo removeing ldradius and coordoutput
    rm ${d}/*ldradius*
    rm ${d}/coord_output
  fi

  if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.ldpred.${i}.ss ]; then
	echo no ldpred score already exists
	present_cord=`ls ${d}/*ldradius* | wc -l`

	if [ $present_cord -eq 0 ];then
		echo there is not already a ldradius file present

		if [ $reftype == "TGP" ];then
		#1000 genomes
		echo went with TGP
			ldpred coord --gf ~/athena/refs/1000genomes/eur.${chr} --ssf temp_files/ss.${low_author}.${chr} --ssf-format CUSTOM --rs RSID --A1 A1 --A2 A2 --eff BETA --pos BP --chr CHR --pval P --se SE --ncol ESS --eff_type LOGOR --out ${d}/coord_output

		else
		#UKBB
		echo went with UKB
			head -5000 geno_files/${low_author}.${chr}.fam | cut -f1 > ${d}/subset_inds
			plink --bfile geno_files/${low_author}.${chr} --keep-fam ${d}/subset_inds --make-bed --out ${d}/for_ldpred
			ldpred coord --gf ${d}/for_ldpred --ssf temp_files/ss.${low_author}.${chr} --ssf-format CUSTOM --rs RSID --A1 A1 --A2 A2 --eff BETA --pos BP --chr CHR --pval P --se SE --ncol ESS --eff_type LOGOR --out ${d}/coord_output
		fi
	fi

	ldpred gibbs --cf ${d}/coord_output --ldf ${d}/ld_file --h2 $h2 --ldr $rad --f $fval --out ${d}/ldpred.${i}.res
        new_file=`ls  ${d}/ldpred.${i}.res_LDpred_p* | wc -l`
        if [ $new_file -gt 0 ];then
	  echo found new file
          new_file_name=`ls  ${d}/ldpred.${i}.res_LDpred_p*`
	  echo it is $new_file_name
          Rscript helper_scripts/ldpred_beta_switch.R dir$dir $new_file_name $author $chr $i
        fi
	prev_reftpye=$reftype

  fi

  let i=i+1
done
