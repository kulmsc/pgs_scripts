goOn=False
while [ $goOn == "False" ]; do
	openSlots=`ls temp_files/rsids* | wc -l`
	if [ $openSlots -lt 8 ]; then
                sizedir=`du temp_files/ | cut -f1`
		if [ $sizedir -lt 50164096 ];then
			echo NOW WE CAN GO - 2
			echo 2 > temp_files/poss_go
		fi
		sleep $(( ( RANDOM % 30 )  + 1 ))
	else
		echo MUST WAIT FOR ROOM TO GO - 2
		sleep $(( ( RANDOM % 30 )  + 1 ))
	fi
done

