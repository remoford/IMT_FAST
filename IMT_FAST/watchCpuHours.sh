#!/bin/bash

while true
do
	./stats.sh
	qstat
	qhost | head -n1
	qhost | sort -n -k7 | less | tail -n$(( $(qstat | wc -l) - 2))
	sleep 10
	clear
done
