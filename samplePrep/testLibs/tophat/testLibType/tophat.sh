#!/bin/bash

#SBATCH --job-name tophat
#SBATCH --cpus-per-task 4
#SBATCH --output /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/OSNN3R3/tophat.out.txt
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --ntasks 1
#SBATCH --ntasks-per-core 1

THEN=`date`

###############
### NOTE!!! ###
###############

### It appears the reads are incorrectly split by the fastq-dump command, so we 
### need to supply them to tophat in reverse order and use --library-type 
### fr-secondstrand, e.g. tophat [options] reads.2.fastq reads.1.fastq. The use
### of fr-secondstrand is based on visual inspection of the mappings in
### /home/tom/Data/zhang_data/subset/tophat, e.g. at LOC_Os02g38920 and
### LOC_Os12g44350. For how these regions were found see
### /home/tom/Data/zhang_data/subset/htseq-count/compareNumberOfCounts.R

OUTDIR=/home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/OSNN3R3
GENOME=/home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/genome/Osativa_204_v7.0
INDEX=/home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/genome/known/known
READS1=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OSNIP_S1_L001_R1_001.trimmed.fastq
READS2=/home/tom/Desktop/5accessions/samplePrep/testLibs/cutadapt/OSNIP_S1_L001_R2_001.trimmed.fastq


srun -c 4 -n 1 -J tophat \
	tophat -p 4 -o $OUTDIR.fs --mate-inner-dist 131 \
	--mate-std-dev 61 --max-intron-length 5000 --library-type fr-firststrand --no-mixed \
	--transcriptome-index $INDEX $GENOME \
	$READS1 $READS2 

srun -c 4 -n 1 -J tophat \
	tophat -p 4 -o $OUTDIR.ss --mate-inner-dist 131 \
	--mate-std-dev 61 --max-intron-length 5000 --library-type fr-secondstrand --no-mixed \
	--transcriptome-index $INDEX $GENOME \
	$READS1 $READS2 

srun -c 4 -n 1 -J tophat \
	tophat -p 4 -o $OUTDIR.fsrev --mate-inner-dist 131 \
	--mate-std-dev 61 --max-intron-length 5000 --library-type fr-firststrand --no-mixed \
	--transcriptome-index $INDEX $GENOME \
	$READS2 $READS1

srun -c 4 -n 1 -J tophat \
	tophat -p 4 -o $OUTDIR.ssrev --mate-inner-dist 131 \
	--mate-std-dev 61 --max-intron-length 5000 --library-type fr-secondstrand --no-mixed \
	--transcriptome-index $INDEX $GENOME \
	$READS2 $READS1

wait

NOW=`date`
MESSAGE=`cat $OUTDIR/tophat.out.txt`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted $THEN finished\n\n### stdout and stderr from job ###\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr

rm $OUTDIR/tophat.out.txt