
maxGo=3
rm temp_files/poss_go
cat ~/athena/doc_score/qc/cv_files/train_eid.0.2.txt | cut -f1 > temp_files/brit_eid
cat ~/athena/doc_score/qc/cv_files/test_eid.0.8.txt | cut -f1 >> temp_files/brit_eid

for (( i=1; i<=$maxGo; i++ )); do
  echo $i >> temp_files/poss_go
done

echo 1 > temp_files/counter

ls refine_other_ss/ | while read scorename;do
    for cchr in {1..22};do

        counter_var=`cat temp_files/counter`
        echo $counter_var
        ./simple_score.sh $scorename $cchr $counter_var &> logs/log.${counter_var}.log &
        echo $counter_var + 1 | bc > temp_files/counter

        sleep $(( ( RANDOM % 10 )  + 1 ))

        goOn=False
        while [ $goOn == "False" ]; do
          openSlots=`cat temp_files/poss_go | wc -l`
          sizedir=`du temp_files/ | cut -f1`
          if [ $openSlots -gt 0 ]; then
            if [ $sizedir -lt 20164096 ];then
              echo NOW WE CAN GO
              goOn=True
            fi
          else
            echo MUST WAIT FOR ROOM TO GO
            sleep $(( ( RANDOM % 5 )  + 1 ))
          fi
        done

  done
done
