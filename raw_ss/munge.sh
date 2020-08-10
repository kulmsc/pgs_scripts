cat common_files/list_authors | while read author;do
  low_author=`echo "$author" | tr '[:upper:]' '[:lower:]'`
  munge_sumstats --sumstats ${author}/clean_${low_author}.txt.gz --N-col ESS --out ${author}/${low_author}.munged
done
