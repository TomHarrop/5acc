#!/bin/bash

for bamfile in $(ls /data/projects/evoreprice/Tom/star/results/*.bam)
do
	library_name=$(basename $bamfile .Aligned.out.bam)
	echo -e "[library_name:]\t$library_name\n[bamfile:]\t$bamfile"
	echo -e "qsub -v bamfile=$bamfile,library_name=$library_name /data/projects/evoreprice/Tom/scripts/htseqScript.sh"
done
