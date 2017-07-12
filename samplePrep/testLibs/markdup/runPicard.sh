#!/bin/bash

#SBATCH --job-name picard
#SBATCH --cpus-per-task 4
#SBATCH --output /home/tom/Desktop/5accessions/samplePrep/testLibs/markdup/picard.out.txt
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --ntasks 1

THEN=`date`

for bamfile in $(ls /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/*/accepted_hits.bam)
do
#	echo $bamfile
	lib=`echo $bamfile | awk -F "/" '{print $9}'`
	lib=`basename $lib .tophat`
#	echo $lib
	#srun -c 4 -J $lib java -jar /usr/local/bin/picard/MarkDuplicates.jar INPUT=$bamfile OUTPUT=$lib.markdup.bam METRICS_FILE=$lib.markdup.metrics 
	srun -c 4 -J $lib java -jar /usr/local/bin/picard/CollectInsertSizeMetrics.jar INPUT=$bamfile OUTPUT=$lib.ISM.txt HISTOGRAM_FILE=$lib.ISM.pdf
done
wait

NOW=`date`
MESSAGE=`cat /home/tom/Desktop/5accessions/samplePrep/testLibs/markdup/picard.out.txt`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted at $THEN\n\nSTDOUT FROM JOB:\n$MESSAGE" | msmtp thomas.harrop@ird.fr

rm /home/tom/Desktop/5accessions/samplePrep/testLibs/markdup/picard.out.txt