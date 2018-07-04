#!/bin/bash

while true
do
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
	du -sh output
	echo "Finished: $COMPLETED_ONESTAGE onestage, $COMPLETED_TWOSTAGE twostage, $COMPLETED_THREESTAGE threestage"
	echo "Points: $POINTS_ONESTAGE onestage, $POINTS_TWOSTAGE twostage, $POINTS_THREESTAGE, threestage"
	qstat -j $(ls output/latest -l | tr '/' '\n' | tail -n1) | grep usage
	#qhost
	sleep $((60*5))
done
