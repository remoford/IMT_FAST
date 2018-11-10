#!/bin/bash
#$ -t 1-18
#$ -e /dev/null
#$ -o /dev/null
task_index=$(( $SGE_TASK_ID - 1  ))
data_sets=(AT1  CHX  DMSO  erlotinib  FUCCI  MCF)
modes=(onestage twostage threestage)
datestamp="_""$(date +%s)"
if [$JOB_NAME -eq ""]
then
	out_dir="output/""$datestamp"
	mkdir $out_dir
	rm output/latest
	ln -s $(pwd)/$out_dir $(pwd)/output/latest
	qsub -N $datestamp -q intel.q -cwd $0 -e /dev/null -o /dev/null #-l hostname=compute-1-[0-19]
# hostname="compute-0-0.local+compute-0-10.local+compute-0-11.local+compute-0-12.local+compute-0-14.local+compute-0-15.local+compute-0-16.local+compute-0-17.local+compute-0-19.local+compute-0-1.local+compute-0-2.local+compute-0-3.local+compute-0-4.local+compute-0-5.local+compute-0-6.local+compute-0-7.local+compute-0-9.local"
	sleep 10
	qstat
	find $out_dir | grep out.txt | xargs -n1 head -n1
else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$(pwd)
	data=${data_sets[ $(( $task_index % 6 )) ]}
	mode=${modes[ $(( $task_index / 6 )) ]}
	binary="./icc_fast"
	output_dir="output/""$JOB_NAME""/"
	(
		echo -n "Job ""$JOB_ID"".""$task_index"" started on ""$HOSTNAME"" from ""$SGE_O_HOST"" logging to ""$output_dir"
		echo
		date
		echo
		uptime
		echo
		sed 's/^$/%/' < /proc/cpuinfo | tr '\n' '^' | tr '%' '\n' | head -n1 | tr '^' '\n'
		echo $(grep "^processor" /proc/cpuinfo | wc -l)" threads"
		echo
		free -h
		echo
		df -h
		echo
		echo "Running $binary $mode $data ..."
		echo
		time ($binary analysis $mode "data/""$data"".txt" | grep -v "range error")
	) 2>&1 > "$output_dir""$binary""_""$mode""_""$data""_out.txt"
fi
