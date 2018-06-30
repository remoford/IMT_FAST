#!/bin/bash

datestamp=$(date +%s)

mkdir "output/""$datestamp"
pushd "output/""$datestamp"

pwd

for mode in onestage twostage threestage
do
	for data in  ../../data/*
	do
		../../icc_fast $mode $data > "icc_""$mode""_""$(basename $data | sed 's/\.txt$//')""_out.txt" &
		../../gcc_fast $mode $data > "gcc_""$mode""_""$(basename $data | sed 's/\.txt$//')""_out.txt" &
	done
done


wait
popd
