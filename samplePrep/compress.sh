#!/bin/bash

#SBATCH --job-name compress
#SBATCH --cpus-per-task 1
#SBATCH --nice=900
#SBATCH --ntasks 1

FILEDIR='/home/tom/Desktop/5accessions/samplePrep'

THEN=`date`
BEFORE=`du -h $FILEDIR --summarize`
find $FILEDIR -name *fastq -type f -exec bash -c 'srun gzip --best {}' \;

# for file in /home/tom/Desktop/laserdissect/*.fastq 
# do
# 	srun gzip --best $file
# done

wait

NOW=`date`
AFTER=`du -h $FILEDIR --summarize`
MESSAGE="GZIP jobs for $FILEDIR finished"
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted $THEN finished\n\n$MESSAGE\n\nDisk usage before:\n\n$BEFORE\n\nDisk usage after:\n\n$AFTER" | msmtp thomas.harrop@ird.fr