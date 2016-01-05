#!/bin/bash

set -e

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

echo -e "[ "$(date)": Genome generation with STAR ]"

# make output directory
outdir="output/star-index"
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# set parameters
genomeFastaFiles="data/genome/os/Osativa_204_v7.0.fa"
sjdbGTFfile="data/genome/os/Osativa_204_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf"
sjdbGTFtagExonParentTranscript="oId"
sjdbGTFtagExonParentGene="gene_name"
sjdbOverhang=109
outFileNamePrefix="$outdir/"

cat << _EOF_
	[ $(date): Submitting job ]
				  genomeFastaFiles:  $genomeFastaFiles
					   sjdbGTFfile:  $sjdbGTFfile
	sjdbGTFtagExonParentTranscript:  $sjdbGTFtagExonParentTranscript
					  sjdbOverhang:  $sjdbOverhang
_EOF_

FAIL=0
cmd="STAR --runThreadN "$maxCpus" --runMode genomeGenerate --genomeDir $outdir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbGTFtagExonParentTranscript $sjdbGTFtagExonParentTranscript --sjdbGTFtagExonParentGene $sjdbGTFtagExonParentGene --sjdbOverhang $sjdbOverhang --outFileNamePrefix $outFileNamePrefix"
srun --output $outdir/stargg.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &

echo -e "[ "$(date)": Waiting for jobs to finish ]"
fail_wait

# log metadata
version="$(STAR --version)"
cat -t <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	STAR version,$version
	genomeFastaFiles,$genomeFastaFiles
	sjdbGTFfile,$sjdbGTFfile
	sjdbOverhang,$sjdbOverhang
_EOF_
echo -e "[ "$(date)": Jobs finished, exiting ]"

exit 0