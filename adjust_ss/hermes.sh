chr=$1
method=$2
author=$3

echo start hermes
echo chr $chr
echo method $method
echo author $author

dir=`head -1 temp_files/poss_dirs`
grep -v -w $dir temp_files/poss_dirs > temp_files/temp_poss
mv temp_files/temp_poss temp_files/poss_dirs

echo dir $dir
let comp2=dir+15
let comp3=dir+30

if [ $method == "clump" ]; then
  echo going clump
  taskset -c ${dir},${comp2},${comp3} ./${method}.sh $chr $author $dir

elif [ $method == "ldpred" ]; then
  echo going ldpred
  taskset -c ${dir},${comp2},${comp3} ./${method}.sh $chr $author $dir

elif [ $method == "ldpred2" ]; then
  echo going ldpred2
  taskset -c ${dir},${comp2},${comp3} ./${method}.sh $chr $author $dir

elif [ $method == "sblup" ]; then
  echo going sblup
  taskset -c ${dir},${comp2},${comp3} ./${method}.sh $chr $author $dir

elif [ $method == "double_weight" ]; then
  echo going sblup
  taskset -c ${dir},${comp2},${comp3} ./${method}.sh $chr $author $dir

fi

rm -r comp_zone/dir${dir}/*
echo yes >> done_check/${author}.${chr}.done
echo $dir >> temp_files/poss_dirs
