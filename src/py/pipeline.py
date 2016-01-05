#!/usr/bin/python3

# slurm options
#SBATCH --ntasks=1
#SBATCH --job-name="pipeline"
#SBATCH --mail-type=ALL
#SBATCH --output=ruffus/pipeline.%j.log

# imports

from ruffus import *
import ruffus.cmdline as cmdline
import datetime
import os
import re
from subprocess import Popen, PIPE

# command-line options
parser = cmdline.get_argparse(description = 'Run LMD analysis pipeline.')
parser.add_argument('--email', '-e',
                        help ='Logon email address for JGI',
                        type = str,
                        dest = 'jgiLogon')
parser.add_argument('--password', '-p',
                        help ='JGI password',
                        type = str,
                        dest = 'jgiPassword')
options = parser.parse_args()
jgiLogon = options.jgiLogon
jgiPassword = options.jgiPassword

# parse SLURM job-id
if 'SLURM_JOBID' in os.environ:
    slurm_jobid = os.environ.get('SLURM_JOBID')
else:
    slurm_jobid = 'local-' + str(datetime.date.today()) 

# time function
def print_now():
    now = datetime.datetime.now()
    return(now.strftime("%Y-%m-%d %H:%M"))

# Custom job submission step. Dirty hack. Need to exit 0 at end of each script
# and use set -e for safety
def submit_job(jobScript, ntasks, cpus_per_task, job_name):
    '''
    Submit the job using salloc hack. When complete return job id and write output to file.
    '''
    # call salloc as subprocess
    proc = Popen(['salloc', '--ntasks=' + ntasks, '--cpus-per-task=' + cpus_per_task,
    '--job-name=' + job_name, jobScript], stdout = PIPE, stderr = PIPE)
    # get stdout and stderr    
    out, err = proc.communicate()
    # parse stderr (salloc output) for job id
    jobRegex = re.compile(b'\d+')
    jobIdBytes = jobRegex.search(err).group(0)
    jobId = jobIdBytes.decode("utf-8")
    # write stderr & stdout to log file    
    outFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.out.txt'
    f = open(outFile, 'wb')
    f.write(out)
    f.close()
    errFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.err.txt'
    f = open(errFile, 'wb')
    f.write(err)
    f.close()
    # mail output
    if proc.returncode != 0:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " FAILED"
    else:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " finished"
    mail = Popen(['mail', '-s', subject, '-A', outFile, '-A', errFile, 'tom'], stdin = PIPE)
    mail.communicate()
    # check subprocess exit code
    os.remove(outFile)    
    os.remove(errFile)
    assert proc.returncode == 0, "Job " + job_name + " failed with non-zero exit code"
    return(jobId)

# Variant of submit_job for JGI jobs where email/password is required
def submit_download_job(jobScript, job_name, jgiLogon, jgiPassword):
    # make sure logon and password were set
    assert jgiLogon, "--email is required"
    assert jgiPassword, "--password is required"

    ntasks = '1'

    # call download script    
    proc = Popen(['salloc', '--ntasks=' + ntasks,'--job-name=' + job_name,
                  jobScript, "-e", jgiLogon, "-p", jgiPassword],
                  stdout = PIPE, stderr = PIPE)
    # get stdout and stderr    
    out, err = proc.communicate()
    # parse stderr (salloc output) for job id
    jobRegex = re.compile(b'\d+')
    jobIdBytes = jobRegex.search(err).group(0)
    jobId = jobIdBytes.decode("utf-8")
    # write stderr & stdout to log file    
    outFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.out.txt'
    f = open(outFile, 'wb')
    f.write(out)
    f.close()
    errFile = 'ruffus/' + job_name + '.' + jobId + '.ruffus.err.txt'
    f = open(errFile, 'wb')
    f.write(err)
    f.close()
    # mail output
    if proc.returncode != 0:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " FAILED"
    else:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " finished"
    mail = Popen(['mail', '-s', subject, '-A', outFile, '-A', errFile, 'tom'], stdin = PIPE)
    mail.communicate()
    # check completion    
    os.remove(outFile)    
    os.remove(errFile)
    assert proc.returncode == 0, "Job " + job_name + " failed with non-zero exit code"
    return(jobId)

# touch function for updating ruffus flag files
def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
            dir_fd=None if os.supports_fd else dir_fd, **kwargs)


# get default pipeline
main_pipeline = Pipeline.pipelines["main"]

#---------------------------------------------------------------
# SETUP TASKS
#---------------------------------------------------------------

#---------------------------------------------------------------
# download rice genome
#

def downloadGenome_sh(outputFiles, jgiLogon, jgiPassword):
    jobScript = 'src/sh/downloadGenome.sh'
    job_name = 'downloadGenome_sh'
    jobId = submit_download_job(jobScript, job_name, jgiLogon, jgiPassword)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    
os_genome = main_pipeline.originate(downloadGenome_sh, "data/genome/os/METADATA.csv", jgiLogon, jgiPassword)

#---------------------------------------------------------------
# define reads
#

def defineOsReads(outputFiles):
    pathToReads = 'data/reads/os'
    assert os.path.isdir(pathToReads), "Error: reads folder " + pathToReads + " missing"
    readFiles = os.listdir(pathToReads)
    print("[", print_now(), ": Using Oryza sativa reads in folder " + pathToReads + " ]")
    for fileName in readFiles:
        qualName = pathToReads + '/' + fileName
        assert os.path.isfile(qualName), "Error: read file " + qualName + " missing"
        print(qualName)
    touch(outputFiles)

os_reads = main_pipeline.originate(defineOsReads, "ruffus/os_reads")

#---------------------------------------------------------------
# JOB STEPS
#---------------------------------------------------------------

#---------------------------------------------------------------
# trim reads
#

def cutadapt_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/cutadapt.sh'
    ntasks = '7'
    cpus_per_task = '1'
    job_name = 'cutadapt'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

os_reads = main_pipeline.transform(task_func = cutadapt_sh,
                                   input = defineOsReads,
                                   filter = suffix("ruffus/os_reads"),
                                   output = "output/cutadapt/METADATA.csv")

#---------------------------------------------------------------
# RUN THE PIPELINE
#---------------------------------------------------------------

# print the flowchart
pipeline_printout_graph("ruffus/flowchart." + slurm_jobid + ".pdf", "pdf")

# run the pipeline (disabled for now)
cmdline.run(options, multithread = 8)