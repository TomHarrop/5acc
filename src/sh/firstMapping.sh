#!/bin/bash

set -e

# catch species code
species=$1

# bash traceback code from https://docwhat.org/tracebacks-in-bash/

_showed_traceback=f

_exit_trap () {
  local _ec="$?"
  if [[ $_ec != 0 && "${_showed_traceback}" != t ]]; then
    traceback 1
  fi
}

_err_trap() {
  local _ec="$?"
  local _cmd="${BASH_COMMAND:-unknown}"
  traceback 1
  _showed_traceback=t
  echo "The command ${_cmd} exited with exit code ${_ec}." 1>&2
}

traceback() {
  # Hide the traceback() call.
  local -i start=$(( ${1:-0} + 1 ))
  local -i end=${#BASH_SOURCE[@]}
  local -i i=0
  local -i j=0

  echo "Traceback (last called is first):" 1>&2
  for ((i=${start}; i < ${end}; i++)); do
    j=$(( $i - 1 ))
    local function="${FUNCNAME[$i]}"
    local file="${BASH_SOURCE[$i]}"
    local line="${BASH_LINENO[$j]}"
    echo "     ${function}() in ${file}:${line}" 1>&2
  done
}

# traps
trap _err_trap SIGHUP SIGINT SIGTERM
trap _exit_trap EXIT
trap _err_trap ERR

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

# how many CPUs we got?
if [[ $SLURM_JOB_CPUS_PER_NODE ]]; then
  maxCpus="$SLURM_JOB_CPUS_PER_NODE"
  echo -e "[ "$(date)": Running with "$maxCpus" CPUs ]"
else
  maxCpus=1
fi

### CODE STARTS HERE ------------------------------------------------------------------

set -u

echo -e "[ "$(date)": First mapping step with STAR ]"

# stop if there is no STAR index
star_index_dir="output/star-index"
if [[ ! -d "$star_index_dir" ]]; then
	echo -e "[ "$(date)": No STAR index found ]"
	exit 1
fi
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# stop if there is no cutadapt folder
cutadapt_dir="output/"$species"/cutadapt"
if [[ ! -d "$cutadapt_dir" ]]; then
	echo -e "[ "$(date)": No cutadapt folder found ]"
	exit 1
fi
echo -e "[ "$(date)": Using cutadapt folder $cutadapt_dir ]"

# make output directory
outdir="output/"$species"/STAR"
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# make directory for step 1
if [[ ! -d $outdir/step1 ]]; then
	mkdir -p $outdir/step1
fi

echo -e "[ "$(date)": Submitting step 1 mapping jobs ]"

# STAR options
OPTIONS="--runThreadN "$maxCpus" --genomeDir "$star_index_dir" --outSJfilterReads Unique --readFilesCommand zcat"

# find the R1.fastq.gz files and match the R2 files to run STAR
shopt -s nullglob
fastq_files=("$cutadapt_dir/*R1.fastq.gz")
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
	cmd="STAR $OPTIONS --readFilesIn $fwd_read_file $rev_read_file --outFileNamePrefix $outdir/step1/$library_name. --genomeLoad LoadAndKeep"
	srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &	
done

echo -e "[ "$(date)": Waiting for step 1 jobs to finish ]"
fail_wait

# remove step 1 bamfiles (waste of space)
echo -e "[ "$(date)": Removing alignment files ]"
rm $outdir/step1/*.bam
rm $outdir/step1/*.sam

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

echo -e "[ "$(date)": Jobs finished, exiting ]"

exit 0