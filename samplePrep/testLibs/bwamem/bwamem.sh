#!/bin/bash

THEN=`date`

GENOME=/home/tom/Desktop/senseTest/phytozome204genome/Osativa_204_v7.0.fa
OUTDIR=/home/tom/Desktop/5accessions/samplePrep/testLibs/bwamem

srun -c 2 -J bwamem -e $OUTDIR/Ob.stderr.txt -o $OUTDIR/ob.bwamem.sam \
	bwa mem -t 2 -w 500	$GENOME	/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OB_S2_L001_R1_001.trimmed.fastq /home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OB_S2_L001_R2_001.trimmed.fastq &

srun -c 2 -J bwamem -e $OUTDIR/osn.stderr.txt -o $OUTDIR/osn.bwamem.sam \
	bwa mem -t 2 -w 500	$GENOME	/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OSNIP_S1_L001_R1_001.trimmed.fastq /home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OSNIP_S1_L001_R2_001.trimmed.fastq &

wait

NOW=`date`

MESSAGE=`cat $OUTDIR/Ob.stderr.txt $OUTDIR/osn.stderr.txt`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted $THEN finished\n\n### stdout and stderr from job ###\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr

