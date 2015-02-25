#!/bin/bash

#SBATCH --job-name htseq
#SBATCH --cpus-per-task 4
#SBATCH --nice=500
#SBATCH --ntasks 1
#SBATCH --ntasks-per-core 1
#SBATCH -o /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/stdout.txt
#SBATCH --open-mode=append


THEN=`date`

GTF='/home/tom/Desktop/laserdissect/output_phytozome_no_rRNA/genome/Osativa_204_v7.0.gene_exons.cuffcomp.rRNA_removed.gtf'
TOPHATDIR='/home/tom/Desktop/5accessions/samplePrep/testLibs/tophat'
OUTDIR='/home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count'

mkdir $OUTDIR/tmp

for bamfile in $(ls $TOPHATDIR/*/accepted_hits.bam);
do
	echo $bamfile
	n=`echo $bamfile | awk -F "/" '{ print $9 }'`
	echo $n
	srun -c 4 -n 1 \
		samtools sort -@ 4 -n -o $OUTDIR/tmp/$n.accepted_hits.sorted.bam -O bam -T $OUTDIR/tmp/$n.accepted_hits.sorted.tmp $bamfile
	srun -c 4 -n 1 -o $OUTDIR/$n.htseq-count \
		-e /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/stdout.txt \
		htseq-count -f bam -r name --stranded=reverse -i gene_name $OUTDIR/tmp/$n.accepted_hits.sorted.bam $GTF
done

# srun -c 4 -n 1 \
# 	samtools sort -@ 4 -n -o /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.fs/accepted_hits.sorted.bam \
# 	-O bam -T /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.fs/accepted_hits.tmp /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.fs/accepted_hits.bam

# srun -c 4 -n 1 \
# 	samtools sort -@ 4 -n -o /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.ss/accepted_hits.sorted.bam \
# 	-O bam -T /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.ss/accepted_hits.tmp /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.ss/accepted_hits.bam

# srun -c 4 -n 1 -o /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/fsrev.htseq-count \
# 	-e /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/htseq-run.err.txt --open-mode=append \
# 	htseq-count -f bam -r name --stranded=reverse -i gene_name /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.fsrev/accepted_hits.sorted.bam $GTF

# srun -c 4 -n 1 -o /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/ssrev.htseq-count \
# 	-e /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/htseq-run.err.txt --open-mode=append \
# 	htseq-count -f bam -r name --stranded=reverse -i gene_name /home/tom/Desktop/5accessions/samplePrep/testLibs/tophat/testLibType/OSNN3R3.ssrev/accepted_hits.sorted.bam $GTF

NOW=`date`
MESSAGE=`cat /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/stdout.txt`
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job finished at $NOW\nJob submitted $THEN finished\n\n### stdout and stderr from job ###\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr

rm /home/tom/Desktop/5accessions/samplePrep/testLibs/htseq-count/stdout.txt

rm -r $OUTDIR/tmp