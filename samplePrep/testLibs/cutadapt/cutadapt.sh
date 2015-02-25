#!/bin/bash

THEN=`date`

READ1=/home/tom/Desktop/5accessions/samplePrep/testLibs/OSNIP_S1_L001_R1_001.fastq
READ2=/home/tom/Desktop/5accessions/samplePrep/testLibs/OSNIP_S1_L001_R2_001.fastq
ADAPTER_FWD="TruSeq_Indexed_Adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_REV="TruSeq_Universal_Adapter_rc=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
OUTDIR=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt
read1_bn=`basename $READ1 .fastq`
read1_tmp=$OUTDIR/$read1_bn.tmp.fastq
read2_bn=`basename $READ2 .fastq`
read2_tmp=$OUTDIR/$read2_bn.tmp.fastq

srun -c 1 -J $read1_bn -o $OUTDIR/$read1_bn.report.txt \
	cutadapt -q 20 -m 50 -a $ADAPTER_FWD -o $read1_tmp -p $read2_tmp $READ1 $READ2

srun -c 1 -J $read2_bn -o $OUTDIR/$read2_bn.report.txt \
	cutadapt -q 20 -m 50 -a $ADAPTER_REV -o $OUTDIR/$read2_bn.trimmed.fastq -p $OUTDIR/$read1_bn.trimmed.fastq $read2_tmp $read1_tmp

rm $read1_tmp
rm $read2_tmp

NOW=`date`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\ncutadapt jobs submitted at $THEN finished." | msmtp thomas.harrop@ird.fr
