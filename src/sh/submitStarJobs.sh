#!/bin/bash

for fwd_read_file in $(ls /data/projects/evoreprice/Data/O*/Trimmed/*-cutadapt_cluster-20150522/*R1*fastq)
do
	n=$(basename $fwd_read_file)
	library_name=${n:0:2}
	rev_read_file=$(ls ${fwd_read_file/$library_name*R1/$library_name*R2})
	echo $library_name
	echo "R1: "
	echo $fwd_read_file
	echo "R2: "
	echo $rev_read_file
	qsub -v fwd_read_file=$fwd_read_file,rev_read_file=$rev_read_file,library_name=$library_name /data/projects/evoreprice/Tom/scripts/starScript.sh
done
