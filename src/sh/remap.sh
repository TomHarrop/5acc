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

set -u

### CODE STARTS HERE ------------------------------------------------------------------

# make output directory
outdir="output/"$species"/STAR/remap"
if [[ ! -d "$outdir" ]]; then
  mkdir -p "$outdir"
fi

echo -e "[ "$(date)": Remapping unmapped reads with STAR ]"

# stop if there is no STAR index
star_index_dir="output/star-index"
if [[ ! -d "$star_index_dir" ]]; then
  echo -e "[ "$(date)": No STAR index found ]"
  exit 1
fi
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# stop if there is no step 1 folder
step1_dir="output/"$species"/STAR/step1"
if [[ ! -d "$step1_dir" ]]; then
  echo -e "[ "$(date)": First step splice junctions not found ]"
  exit 1
fi

# recover the splice junction files from step 1
echo -e "[ "$(date)": Finding splice junctions in "$step1_dir" ]"

shopt -s nullglob
first_pass_junctions=(""$step1_dir"/*SJ.out.tab")
shopt -u nullglob

cat <<- _EOF_
  [ $(date): Found SJ.out.tab files from step 1 ]
  $(for tab in $first_pass_junctions; do echo $tab; done)
_EOF_

# set relaxed mapping options
# setting --outFilterMismatchNmax to a high number means that number of mismatches will be 
# controlled by the --outFilterMismatchNoverLmax parameter (defaults to 0.3, so about 60
# mismatches for a 200nt alignment). Setting --outFilterMatchNminOverLread to 0.45 means
# that 50nt alignments are acceptable for a read of 112nt.

relaxedMapping="--outFilterMismatchNmax 1000 --outFilterMatchNminOverLread 0.45"

# set STAR options
OPTIONS=""$relaxedMapping" --sjdbFileChrStartEnd "$first_pass_junctions" --genomeLoad NoSharedMemory --runThreadN "$maxCpus" --genomeDir "$star_index_dir" --outSJfilterReads Unique --readFilesCommand zcat --outSAMtype BAM Unsorted --quantMode GeneCounts --outBAMcompression 10 --outReadsUnmapped Fastx"

# find unmapped mate1 files
star_dir="output/"$species"/STAR"
shopt -s nullglob
unmappedMate1Files=(""$star_dir"/*.Unmapped.out.mate1.gz")
shopt -u nullglob

# match mate1 to mate2 and map
for fwdMate in $unmappedMate1Files
do
  n=$(basename $fwdMate)
  library_name=${n:0:2}
  revMate=${fwdMate/mate1/mate2}
  # double check rev_read_file exists
  if [[ ! -e $revMate ]]; then
    echo -e "[ "$(date)" : Couldn't find reverse read file $revMate for library $library_name ]"
    exit 1
  fi
  cat <<- _EOF_
  [ $(date) : Submitting STAR run ]
  library_name: $library_name
  fwdMate: $fwdMate
  revMate: $revMate
_EOF_
  cmd="STAR $OPTIONS --readFilesIn $fwdMate $revMate --outFileNamePrefix $outdir/$library_name."
  srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &
done

echo -e "[ "$(date)": Waiting for remapping jobs to finish ]"
fail_wait

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
  Script,${0}
  branch,$(git rev-parse --abbrev-ref HEAD)
  hash,$(git rev-parse HEAD)
  date,$(date +%F)
  STAR version,$(STAR --version)
  STAR index,$star_index_dir
_EOF_

echo -e "[ "$(date)": Jobs finished, exiting ]"

exit 0