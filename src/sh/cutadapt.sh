#!/bin/bash

#SBATCH --job-name cutadapt
#SBATCH --ntasks 6
#SBATCH --output /tmp/cutadapt.%N.%j.out
#SBATCH --open-mode=append

# prep
THEN="$(date)"

# make today's output directory
outdir="output/cutadapt-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
version="$(cutadapt --version)"
cat -t <<- _EOF_ > $outdir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cutadapt version,$version
_EOF_

# catch SIGTERM etc.
clean_up() {
	# remove temp files before exit
	echo -e "[ "$(date)" : Script aborted ]"
	# email output
	NOW="$(date)"
	MESSAGE="$(cat /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
	cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" aborted" tom
	Job $SLURM_JOBID submitted at $THEN was aborted.
	Concatenated stdout files:
	$MESSAGE
_EOF_
	mv /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# parameters
adaptorFwd='TruSeq_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTC'
adaptorRev='Illumina_single_end=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG'
trim_qualities=20
minimum_length=50

echo -e "[ "$(date)": Submitting cutadapt jobs ]"
for fwd_reads in $(ls data/reads/*R1.fastq.gz); do
	fFile="$(basename $fwd_reads)"
	lib_name="${fFile:0:2}"
	rev_reads="data/reads/$(basename $fwd_reads 1.fastq.gz)2.fastq.gz"
	# check that rev_reads are really there
	if [[ ! -e $rev_reads ]]; then
		echo "Error: rev_reads not found\n[ lib_name ]:\t$lib_name\n[ rev_reads ]:\t$rev_reads"
		exit 1
	fi
	output="$outdir/$lib_name.R1.fastq.gz"
	paired_output="$outdir/$lib_name.R2.fastq.gz"
	# print some info
	echo -e "Running cutadapt:\n[ lib_name ]:\t$lib_name\n[ fwd_reads ]:\t$fwd_reads\n[ rev_reads ]:\t$rev_reads\n[ R1 out ]:\t$output\n[ R2 out ]:\t$paired_output"
	# run cutadapt
	cmd="cutadapt -a $adaptorFwd -A $adaptorRev --quality-cutoff=$trim_qualities --minimum-length $minimum_length --output=$output --paired-output=$paired_output $fwd_reads $rev_reads"
	srun --output $outdir/$lib_name.out --exclusive --ntasks=1 --cpus-per-task=1 $cmd &
done

echo -e "[ "$(date)": Waiting for jobs to finish ]"

wait

echo -e "[ "$(date)": Jobs finished, tidying up ]"

# email output
NOW="$(date)"
MESSAGE="$(cat /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID submitted at $THEN is finished.
	Concatenated stdout files:
	$MESSAGE
_EOF_

mv /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/

echo -e "[ "$(date)": Exiting ]"
exit 0
