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

# stop if there is no STAR output
star_dir="output/"$species"/STAR"
if [[ ! -d "$star_dir" ]]; then
	echo -e "[ "$(date)": STAR output not found ]"
	exit 1
fi

# get some stats
cat <<- _EOF_ > "$star_dir"/compressStats.txt
	$(gzip --version)
	[ $(date) : Disk usage before compression ]
	$(du -shc $star_dir/*.Unmapped.out.mate*)
_EOF_

find "$star_dir" -type f -name "*.Unmapped.out.mate*" -exec bash -c \
	'echo -e "[ "$(date)" : compressing {} ]" ;
	srun --exclusive --ntasks=1 --cpus-per-task=1 gzip --best {} &' \;
echo -e "[ "$(date)" : waiting for gzip jobs to finish ]"
fail_wait

# more stats
cat <<- _EOF_ >> "$star_dir"/compressStats.txt
	[ $(date) : Disk usage after compression ]
	$(du -shc $star_dir/*.Unmapped.out.mate*.gz)
_EOF_

exit 0