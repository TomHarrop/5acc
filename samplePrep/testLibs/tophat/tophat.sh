#!/bin/bash

#SBATCH --job-name tophat
#SBATCH --cpus-per-task 5
#SBATCH --output /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/tophat.out.txt
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --ntasks 1
#SBATCH --ntasks-per-core 1

THEN=`date`

GENOME=/home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/genome/Osativa_204_v7.0
INDEX=/home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/genome/known/known
READS1=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OSNIP_S1_L001_R1_001.trimmed.fastq
READS2=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OSNIP_S1_L001_R2_001.trimmed.fastq
READS3=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OB_S2_L001_R1_001.trimmed.fastq
READS4=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OB_S2_L001_R2_001.trimmed.fastq

srun -c 5 -n 1 -J tophat \
	tophat -p 4 --mate-inner-dist 131 \
	--mate-std-dev 61 --max-intron-length 5000 \
	--library-type fr-firststrand --no-mixed \
	-o /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/OSNN3R3 \
	--transcriptome-index $INDEX $GENOME \
	$READS1 $READS2

srun -c 5 -n 1 -J tophat \
	tophat -p 4 --mate-inner-dist 131 \
	--mate-std-dev 61 --max-intron-length 5000 \
	--library-type fr-firststrand --no-mixed \
	-o /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/OBN1R2 \
	--transcriptome-index $INDEX $GENOME \
	$READS3 $READS4

wait

NOW=`date`
MESSAGE=`cat /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/tophat.out.txt`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted $THEN finished\n\n### stdout and stderr from job ###\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr

rm /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/tophat.out.txt