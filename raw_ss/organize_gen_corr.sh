rm genetic_correlations

cat common_files/list_authors | while read fauth1;do
  cat common_files/list_authors | while read fauth2;do

  auth1=`echo "$fauth1" | tr '[:upper:]' '[:lower:]'`
  auth2=`echo "$fauth2" | tr '[:upper:]' '[:lower:]'`

  if [ -e gen_corr/ldsc/${auth1}.${auth2}.corr.log ];then
    rg_ldsc=`cat gen_corr/ldsc/${auth1}.${auth2}.corr.log | fgrep "Genetic Correlation:" | cut -f3 -d' '`
    se_ldsc=`cat gen_corr/ldsc/${auth1}.${auth2}.corr.log | fgrep "Genetic Correlation:" | cut -f4 -d' ' | cut -f2 -d'(' | cut -f1 -d')'`
  else
    rg_ldsc=NA
    se_ldsc=NA
  fi

  if [ -e gen_corr/hdl/${auth1}.${auth2}.corr.log ];then
    rg_hdl=`cat gen_corr/hdl/${auth1}.${auth2}.corr.log | fgrep Correlation | cut -f2 -d' '`
    se_hdl=`cat gen_corr/hdl/${auth1}.${auth2}.corr.log | fgrep Correlation | cut -f3 -d' '`
  else
    rg_hdl=NA
    se_hdl=NA
  fi

  echo $auth1 $auth2 $rg_ldsc $se_ldsc $rg_hdl $se_hdl >> genetic_correlations

  done
done
