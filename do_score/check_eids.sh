echo hello > all_eid
ls small_score_files/*zst | while read line;do
  zstd -dc ${line} | sed 's/  */\t/g' | cut -f2 > temp
  paste all_eid temp > temp2; mv temp2 all_eid
done
