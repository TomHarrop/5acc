
#!/bin/bash

for bamfile in $(ssh ***REMOVED*** "ls /data/projects/evoreprice/Tom/star/results/*.bam")
do
	echo $bamfile
	library_name="$(basename $bamfile .Aligned.out.bam)"
	echo -e "/home/tom/Desktop/5accessions/src/sh/runHtseqLocal.sh $bamfile $library_name &"
done