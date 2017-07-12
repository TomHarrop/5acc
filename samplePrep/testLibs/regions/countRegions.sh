#!/bin/bash

#SBATCH --job-name regions2
#SBATCH --cpus-per-task 1
#SBATCH --output /home/tom/Desktop/5accessions/samplePrep/testLibs/regions/Ob/countRegions.out.txt
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --ntasks 1

OUTDIR='/home/tom/Desktop/5accessions/samplePrep/testLibs/regions/OsNIP'
INFILE='/home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/OSNN3R3/accepted_hits.bam'
SPECIES='OsNipn3r3'

THEN=`date`

## split rRNA and non-rRNA reads into separate files

srun -c 1 -n 1 \
	samtools view -b -o $OUTDIR/$SPECIES.rRNA.bam -U $OUTDIR/$SPECIES.norRNA.bam -b -L /home/tom/Data/rRNA/rRNAs2.named.noheader.bed \
	$INFILE

## split non-rRNA reads into reads overlapping the GTF and reads not overlapping the GTF

srun -c 1 -n 1 \
	samtools view -b -o $OUTDIR/$SPECIES.norRNA.inGTF.bam -U $OUTDIR/$SPECIES.norRNA.notInGTF.bam -b \
	-L /home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/RSeqC/Osativa_204_v7.gene_exons.cuffcomp.rRNA_removed.bed \
	$OUTDIR/$SPECIES.norRNA.bam

## reads mapped in pairs (primary alignments only) = total proper pairs

srun -c 1 -n 1 -o $OUTDIR/properPairs.txt \
	samtools view -c -f 67 -F 256 $INFILE

## total proper paired pairs that came from rRNA regions

srun -c 1 -n 1 -o $OUTDIR/rRNApairs.txt \
	samtools view -c -f 67 -F 256 $OUTDIR/$SPECIES.rRNA.bam

## total proper pairs that are not rRNA and are not in genes ('intergenic' / ncRNAs?)

srun -c 1 -n 1 -o $OUTDIR/notInExons.txt \
	samtools view -c -f 67 -F 256 $OUTDIR/$SPECIES.norRNA.notInGTF.bam


## total proper pairs that are not rRNA and WERE in the GTF

srun -c 1 -n 1 -o $OUTDIR/inExons.txt \
	samtools view -c -f 67 -F 256 $OUTDIR/$SPECIES.norRNA.inGTF.bam

## try to make a nice tab-delimited table for printing to --output file

properPairs=`cat $OUTDIR/properPairs.txt`
rRNApairs=`cat $OUTDIR/rRNApairs.txt`
notInExons=`cat $OUTDIR/notInExons.txt`
inExons=`cat $OUTDIR/inExons.txt`

HEADER='properPairs\trRNApairs\tnotInExons\tinExons\n'

printf "$HEADER\t$properPairs\t$rRNApairs\t$notInExons\t$inExons\n" > $OUTDIR/summary.tab

NOW=`date`
MESSAGE=`cat /home/tom/Desktop/5accessions/samplePrep/testLibs/regions/Ob/countRegions.out.txt`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted $THEN finished\n\n### stdout and stderr from job ###\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr

rm /home/tom/Desktop/5accessions/samplePrep/testLibs/regions/Ob/countRegions.out.txt