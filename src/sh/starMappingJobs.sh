#!/bin/bash

#SBATCH --job-name star
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --output /tmp/star.%N.%j.out
#SBATCH --open-mode=append
#SBATCH --nice=500
#SBATCH --mail-user=thomas.harrop@ird.fr
#SBATCH --mail-type=ALL

THEN="$(date)"
echo -e "[ "$(date)": Mapping with STAR ]"

# make today's output directory
outdir="output/STAR-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# choose the most recent STAR index
shopt -s nullglob
folders=(output/star-index*/)
shopt -u nullglob
# stop if there is no STAR index
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": No STAR index found ]"
	exit 1
fi
star_index_dir="${folders[-1]}"
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# choose the most recent cutadapt output
shopt -s nullglob
folders=(output/cutadapt*/)
shopt -u nullglob
# stop if there is no cutadapt folder
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": No cutadapt folder found ]"
	exit 1
fi
cutadapt_dir="${folders[-1]}"
echo -e "[ "$(date)": Using cutadapt folder $cutadapt_dir ]"

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	STAR version,$(STAR --version)
	STAR index,$star_index_dir
	cutadapt folder,$cutadapt_dir
_EOF_

# catch SIGTERM etc
clean_up() {
	echo -e "[ "$(date)" : Script aborted ]"
	# email output
	cat <<- _EOF_ | msmtp thomas.harrop@ird.fr
	To: thomas.harrop@ird.fr
	From: schinkendosen@gmail.com
	Subject: [Tom@SLURM] Job $SLURM_JOBID aborted
	Job $SLURM_JOBID submitted at $THEN was aborted.
	
	Concatenated stdout files:

	$(cat /tmp/star.$SLURM_JOB_NODELIST.$SLURM_JOBID.out)
_EOF_
	mv /tmp/star.$SLURM_JOB_NODELIST.$SLURM_JOBID.out "$outdir"/
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# load genome
echo -e "[ "$(date)": Loading genome into shared memory ]"
cmd="STAR --runThreadN 6 --genomeDir $star_index_dir --genomeLoad LoadAndExit --outFileNamePrefix $outdir/gLoad."
srun --ntasks=1 --exclusive --cpus-per-task=6 $cmd

# find the R1.fastq.gz files and match the R2 files to run STAR
shopt -s nullglob
fastq_files=("$cutadapt_dir*R1.fastq.gz")
shopt -u nullglob
for fwd_read_file in $fastq_files
do
	n=$(basename $fwd_read_file)
	library_name=${n:0:2}
	rev_read_file=${fwd_read_file/$library_name.R1/$library_name.R2}
	# double check rev_read_file exists
	if [[ ! -e $rev_read_file ]]; then
		echo -e "[ "$(date)" : Error! Couldn't find reverse read file $rev_read_file for library $library_name ]"
		exit 1
	fi
	cat <<- _EOF_
	[ $(date) : Submitting STAR run ]
	library_name: $library_name
	fwd_read_file: $fwd_read_file
	rev_read_file: $rev_read_file
_EOF_
	cmd="STAR --runThreadN 6 --genomeDir $star_index_dir --readFilesIn $fwd_read_file $rev_read_file --outFileNamePrefix $outdir/$library_name. --outSAMtype BAM Unsorted --quantMode GeneCounts --genomeLoad LoadAndKeep --readFilesCommand zcat"
	srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task=6 $cmd &	
done

echo -e "[ "$(date)": Waiting for jobs to finish ]"

wait

echo -e "[ "$(date)": Jobs finished, removing index from memory ]"
srun --exclusive --ntasks=1 --cpus-per-task=6 \
	STAR --runThreadN 6 --genomeDir $star_index_dir --genomeLoad Remove --outFileNamePrefix $outdir/gRem.

echo -e "[ "$(date)": Tidying up ]"

# email output
cat <<- _EOF_ | msmtp thomas.harrop@ird.fr
	To: thomas.harrop@ird.fr
	From: schinkendosen@gmail.com
	Subject: [Tom@SLURM] Job $SLURM_JOBID finished
	Job $SLURM_JOBID submitted at $THEN is finished.

	Job log:
	$(cat /tmp/star.$SLURM_JOB_NODELIST.$SLURM_JOBID.out)
_EOF_

mv /tmp/star."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/

echo -e "[ "$(date)": Exiting ]"
exit 0