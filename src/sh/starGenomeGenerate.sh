#!/bin/bash

#SBATCH --job-name stargg
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=6
#SBATCH --output /tmp/stargg.%N.%j.out
#SBATCH --open-mode=append

THEN="$(date)"
echo -e "[ "$(date)": Genome generation with STAR ]"

# make output directory
outdir="output/star-index-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# catch sigkill
clean_up() {
	echo -e "[ "$(date)" : Script aborted ]"
	# email output
	NOW="$(date)"
	MESSAGE="$(cat /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
	printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job "$SLURM_JOBID" aborted\nJob "$SLURM_JOBID" submitted at $THEN was aborted.\n\nConcatenated stdout files:\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr
	mv /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# set parameters
genomeFastaFiles="data/genome/Osativa_204_v7.0.fa"
sjdbGTFfile="data/genome/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"
sjdbGTFtagExonParentTranscript="gene_name"
sjdbOverhang=109

echo -e "[ "$(date)": Submitting job ]\ngenomeFastaFiles:\t\t$genomeFastaFiles\nsjdbGTFfile:\t\t\t$sjdbGTFfile\nsjdbGTFtagExonParentTranscript:\t$sjdbGTFtagExonParentTranscript\nsjdbOverhang:\t\t\t$sjdbOverhang"

cmd="STAR --runThreadN 6 --runMode genomeGenerate --genomeDir $outdir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbGTFtagExonParentTranscript $sjdbGTFtagExonParentTranscript --sjdbOverhang $sjdbOverhang"

srun --output $outdir/stargg.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &

echo -e "[ "$(date)": Waiting for jobs to finish ]"
wait
echo -e "[ "$(date)": Jobs finished, tidying up ]"

# log metadata
version="$(STAR --version)"
echo -e \
"Script\t${0}
STAR version\t$version
date\t$(date +%F)
genomeFastaFiles\t$genomeFastaFiles
sjdbGTFfile\t$sjdbGTFfile
sjdbGTFtagExonParentTranscript\t$sjdbGTFtagExonParentTranscript
sjdbOverhang\t$sjdbOverhang" > $outdir/METADATA.tsv

# email output
NOW="$(date)"
MESSAGE="$(cat /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
printf "To: thomas.harrop@ird.fr\nFrom: schinkendosen@gmail.com\nSubject: [Tom@SLURM] Job "$SLURM_JOBID" finished\nJob "$SLURM_JOBID" submitted at $THEN is finished.\n\nConcatenated stdout files:\n\n$MESSAGE" | msmtp thomas.harrop@ird.fr

mv /tmp/stargg."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/

echo -e "[ "$(date)": Exiting ]"
exit 0