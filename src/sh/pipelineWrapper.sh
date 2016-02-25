#!/bin/bash

while [ "$1" != "" ]; do
	case $1 in
		-e )	shift
				jgi_logon=$1
				;;
		-p )	shift
				jgi_password=$1
				;;
		* )		echo "Bad input"
				exit 1
	esac
	shift
done

set -eu

mail_output() {
	subject="[Tom@SLURM] Pipeline started at $(date) finished"
	echo "" | mail -s "$subject" -A ruffus/pipeline.log.txt tom
	rm ruffus/pipeline.log.txt
}

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
  mail_output
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

# start code
if [[ "$jgi_logon" && "$jgi_password" ]]; then
	src/py/pipeline.py -e "$jgi_logon" -p "$jgi_password" -v5 \
		&> ruffus/pipeline.log.txt
else
	src/py/pipeline.py -v5 &> ruffus/pipeline.log.txt
fi
mail_output
exit 0