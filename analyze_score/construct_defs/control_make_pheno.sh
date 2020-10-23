author=xie
cat parallel_splits | while read line;do
  star=`echo $line | cut -d' ' -f1`
  end=`echo $line | cut -d' ' -f2`
  echo $star
  if [ ! -e raw_output/diag.time.${author}.${star}.txt.gz ];then
    python better_pheno.py $star $end &
    sleep 80
  fi
done
