#!/bin/bash

set -eu

# how many CPUs we got?
if [[ $SLURM_JOB_CPUS_PER_NODE ]]; then
	maxCpus="$SLURM_JOB_CPUS_PER_NODE"
	echo -e "[ "$(date)": Running with "$maxCpus" CPUs ]"
else
	maxCpus=1
fi

# cleanup functions
exit_error() {
	echo -e "[ "$(date)": Script aborted ]"
	exit 1
}

# catch exit codes
trap_exit() {
	exitCode=$?
	if (( "exitCode" == 0 )) ; then
		exit 0
	else
		exit_error
	fi
}

# traps
trap exit_error SIGHUP SIGINT SIGTERM
trap trap_exit EXIT

# handle waiting
FAIL=0
fail_wait() {
for job in $(jobs -p); do
  wait $job || let "FAIL+=1"
done
if [[ ! "$FAIL" == 0 ]]; then
  exit 1
fi
}

### CODE STARTS HERE ------------------------------------------------------------------

# make output directory
outdir="output/cutadapt"
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat -t <<- _EOF_ > $outdir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cutadapt version,$(cutadapt --version)
_EOF_

echo -e "[ "$(date)": Adaptor trimming with cutadapt ]"

# parameters
adaptorFwd='TruSeq_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTC'
adaptorRev='Illumina_single_end=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG'
trim_qualities=20
minimum_length=50

echo -e "[ "$(date)": Submitting cutadapt jobs ]"
fwdReadFiles=("data/reads/os/*R1.fastq.gz")
for fwdReadFile in $fwdReadFiles; do
	fFile="$(basename $fwdReadFile)"
	lib_name="${fFile:0:2}"
	rev_reads="data/reads/os/$(basename $fwdReadFile 1.fastq.gz)2.fastq.gz"
	# check that rev_reads are really there
	if [[ ! -e $rev_reads ]]; then
		echo "Error: rev_reads not found\n[ lib_name ]:\t$lib_name\n[ rev_reads ]:\t$rev_reads"
		exit 1
	fi
	output="$outdir/$lib_name.R1.fastq.gz"
	paired_output="$outdir/$lib_name.R2.fastq.gz"
	echo $output
	echo $paired_output
	# print some info
	echo -e "Running cutadapt:\n[ lib_name ]:\t$lib_name\n[ fwd_reads ]:\t$fwdReadFile\n[ rev_reads ]:\t$rev_reads\n[ R1 out ]:\t$output\n[ R2 out ]:\t$paired_output"
	# run cutadapt
	cmd="cutadapt -a $adaptorFwd -A $adaptorRev --quality-cutoff=$trim_qualities --minimum-length=$minimum_length --output=$output --paired-output=$paired_output $fwdReadFile $rev_reads"
	srun --output $outdir/$lib_name.out --exclusive --ntasks=1 --cpus-per-task=1 $cmd &
done

echo -e "[ "$(date)": Waiting for jobs to finish ]"
fail_wait

echo -e "[ "$(date)": Jobs finished, exiting ]"

exit 0