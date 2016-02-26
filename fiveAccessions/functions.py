#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess
import re
import os
import datetime

############################
# JOB SUBMISSION FUNCTIONS #
############################


def submit_job(job_script, ntasks, cpus_per_task, job_name, extras=[]):
    # type: (str, str, str, str, str, list) -> str
    '''
    Submit the job using salloc hack. When complete return job id and write
    output to file.
    '''
    # call salloc as subprocess
    proc = subprocess.Popen(['salloc', '--ntasks=' + ntasks,
                             '--cpus-per-task=' + cpus_per_task,
                             '--job-name=' + job_name, job_script] + extras,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    # get stdout and stderr
    out, err = proc.communicate()
    # parse stderr (salloc output) for job id
    job_regex = re.compile(b'\d+')
    job_id_bytes = job_regex.search(err).group(0)
    job_id = job_id_bytes.decode("utf-8")
    # write stderr & stdout to log file
    out_file = 'ruffus/' + job_name + '.' + job_id + '.ruffus.out.txt'
    with open(out_file, 'wb') as f:
        f.write(out)
    err_file = 'ruffus/' + job_name + '.' + job_id + '.ruffus.err.txt'
    with open(err_file, 'wb') as f:
        f.write(err)
    # mail output
    if proc.returncode != 0:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " FAILED"
    else:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " finished"
    mail = subprocess.Popen(['mail', '-s', subject, '-A', out_file, '-A',
                             err_file, 'tom'], stdin=subprocess.PIPE)
    mail.communicate()
    # check subprocess exit code
    os.remove(out_file)
    os.remove(err_file)
    assert proc.returncode == 0, ("Job " + job_name +
                                  " failed with non-zero exit code")
    return(job_id)


def print_job_submission(job_name, job_id):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    print('[', now, '] : Job ' + job_name + ' run with JobID ' + job_id)


# touch function for updating ruffus flag files
def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
                 dir_fd=None if os.supports_fd else dir_fd, **kwargs)
