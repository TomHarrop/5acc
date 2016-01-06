#!/bin/bash

set -e

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


# how many CPUs we got?
if [[ $SLURM_JOB_CPUS_PER_NODE ]]; then
  maxCpus="$SLURM_JOB_CPUS_PER_NODE"
  echo -e "[ "$(date)": Running with "$maxCpus" CPUs ]"
else
  maxCpus=1
fi

### CODE STARTS HERE ------------------------------------------------------------------

set -u

echo -e "[ "$(date)": Loading genome into shared memory ]"

# stop if there is no STAR index
star_index_dir="output/star-index"
if [[ ! -d "$star_index_dir" ]]; then
  echo -e "[ "$(date)": No STAR index found ]"
  exit 1
fi
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# make today's output directory
outdir="output/star-index/genomeLoad"
if [[ ! -d $outdir ]]; then
  mkdir -p $outdir
fi

# load genome into memory
cmd="STAR --runThreadN "$maxCpus" --genomeDir $star_index_dir --genomeLoad LoadAndExit --outFileNamePrefix $outdir/gLoad."
srun --ntasks=1 --exclusive --cpus-per-task="$maxCpus" $cmd

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
  Script,${0}
  branch,$(git rev-parse --abbrev-ref HEAD)
  hash,$(git rev-parse HEAD)
  date,$(date +%F)
  STAR version,$(STAR --version)
  STAR index,$star_index_dir
_EOF_

exit 0