chr=$1
author=$2
dir=$3
low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
d=comp_zone/dir${dir}


i=1
cat all_specs/smtpred_param_specs | tail -n +2 | while read spec;do
  max_corr=`echo $spec | cut -f1`
  est_type=`echo $spec | cut -f2`

  if [ ! -e ~/athena/doc_score/mod_sets/${author}/${low_author}.smtpred.${i}.ss ]; then

    if [ $est_type == "normal" ];then
      Rscript helper_scripts/set_up_smtpred.R $low_author $d $max_corr
   
      cat ${d}/ss_files | cut -f1 -d'.' | while read new_author;do
        if [ $new_author == "imsgc" ];then 
          upper_author=`echo $new_author | tr [a-z] [A-Z]`
        else
          upper_author=${new_author^}
        fi
        zcat ~/athena/doc_score/raw_ss/${upper_author}/chr_ss/${new_author}_${chr}.ss.gz > ${d}/${new_author}.ss
      done

      ss_line=`cat ${d}/ss_files | while read line;do echo ${d}/${line}; done | tr '\n' ' '`

      python ~/Programs/smtpred/smtpred.py --betafiles $ss_line --nfile ${d}/samp_size --h2file ${d}/herit --rgfile ${d}/rg

      Rscript helper_scripts/finsih_smtpred.R $d $i $low_author $chr

    fi

    if [ $est_type == "sblup" ];then


    fi

  fi

  let i=i+1
done
