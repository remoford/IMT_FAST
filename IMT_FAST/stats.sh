#!/bin/bash

        HOURS_CPUTIME=$(( $( (grep -E "^ll=" output/latest/* |\
                 tr '\n' '%' |\
                 tr ' ' '\n' |\
                 grep "s$" |\
                 tr -d 's' |\
                 tr '\n' '+' ; echo '0')|\
                         bc |\
                         cut -d'.' -f1 ) / 3600 ))

        COMPLETED_POINTS=$(grep -E "^ll=" output/latest/* | wc -l)

        COMPLETED_ONESTAGE=$(grep "total runtime" output/latest/*onestage* | wc -l)
        COMPLETED_TWOSTAGE=$(grep "total runtime" output/latest/*twostage* | wc -l)
        COMPLETED_THREESTAGE=$(grep "total runtime" output/latest/*threestage* | wc -l)

        POINTS_ONESTAGE=$(grep "^ll=" output/latest/*onestage* | wc -l)
        POINTS_TWOSTAGE=$(grep "^ll=" output/latest/*twostage* | wc -l)
        POINTS_THREESTAGE=$(grep "^ll=" output/latest/*threestage* | wc -l)

        date
        echo -n "$HOURS_CPUTIME hours sum walltime, $COMPLETED_POINTS points "
        pushd output/latest 2>&1 >/dev/null
	du -sh
	popd 2>&1 >/dev/null
        echo "Finished: $COMPLETED_ONESTAGE onestage, $COMPLETED_TWOSTAGE twostage, $COMPLETED_THREESTAGE threestage"
        echo "Points: $POINTS_ONESTAGE onestage, $POINTS_TWOSTAGE twostage, $POINTS_THREESTAGE, threestage"


	for stage in "two" "three"
	do
	echo "$stage"
	grep "finished seedIdx" "output/latest/gcc_fast_""$stage""stage_"* |\
		sed 's#output*stage_##' |\
		sed 's#_out.txt##' |\
		cut -d':' -f1 |\
		cut -d'_' -f4 |\
		sort |\
		uniq -c |\
		sed 's/\ *//' |\
		tr '\n' ' '
		#sed 's/$/\n/'
	echo
	#for file in output/latest/*${stage}stage* ; do echo $file | cut -d'_' -f4 ; echo -n " " ; (cat $file | tr '\n' ' ' | sed 's/finished seedIdx/\n/g' | while read line ; do echo $line | tr ' ' '\n' | grep "iter=" | tail -n1; done | cut -d'=' -f2 | tr '\n' '+' ; echo 0 ) | bc | tr -d '\n' ; echo -ne '   ' ; done | tr -d '\n' ; echo 
	echo
	done
