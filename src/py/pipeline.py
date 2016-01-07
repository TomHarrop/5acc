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
def submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = ""):
    '''
    Submit the job using salloc hack. When complete return job id and write output to file.
    '''
    # call salloc as subprocess
    proc = Popen(['salloc', '--ntasks=' + ntasks, '--cpus-per-task=' + cpus_per_task,
    '--job-name=' + job_name, jobScript, extras], stdout = PIPE, stderr = PIPE)
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
# generate STAR index for OS
#

def starGenomeGenerate_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/starGenomeGenerate.sh'
    ntasks = '1'
    cpus_per_task = '4'
    job_name = 'stargg'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    # update ruffus flag
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

os_index = main_pipeline.transform(task_func = starGenomeGenerate_sh,
                                   input = downloadGenome_sh,
                                   filter = suffix("data/genome/os/METADATA.csv"),
                                   output = "output/star-index/METADATA.csv")


#---------------------------------------------------------------
# MAPPING AND TRIMMING SUBPIPELINE
#---------------------------------------------------------------

#---------------------------------------------------------------
# define reads
#
def defineReads(outputFiles, pathToReads, species):
    assert os.path.isdir(pathToReads), "1 " +outputFiles+ " 2 " + pathToReads + " 3 " + species + "Error: reads folder " + pathToReads + " missing"
    readFiles = os.listdir(pathToReads)
    print("[", print_now(), ": Using reads in folder " + pathToReads + " ]")
    for fileName in readFiles:
        qualName = pathToReads + '/' + fileName
        assert os.path.isfile(qualName), "Error: read file " + qualName + " missing"
        print(qualName)
    touch(outputFiles)

#---------------------------------------------------------------
# trim reads
#
def cutadapt_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/cutadapt.sh'
    ntasks = '2'
    cpus_per_task = '1'
    job_name = 'cutadapt_sh'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = species)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
    
#---------------------------------------------------------------
# load genome into memory
#
def loadGenome_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/loadGenome.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'loadGenome_sh'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# first mapping step (shared memory across jobs)
#
def firstMapping_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/firstMapping.sh'
    ntasks = '1'
    cpus_per_task = '4'
    job_name = 'firstMapping_sh'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = species)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# unload genome from memory
#
def unloadGenome_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/unloadGenome.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'unloadGenome_sh'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# second mapping step (shared memory within job)
#


#---------------------------------------------------------------
# construct the trimming pipeline
#
def makeTrimmingPipeline(pipelineName, pathToReads, species):    
    newPipeline = Pipeline(pipelineName)
    headTask = newPipeline.originate(name = species + " reads",
                                     task_func = defineReads,
                                     output = "ruffus/" + species + ".reads",
                                     extras = [pathToReads, species])
    tailTask = newPipeline.transform(task_func = cutadapt_sh,
                          input = defineReads,
                          filter = suffix("ruffus/" + species + ".reads"),
                          output = "output/" + species + "/cutadapt/METADATA.csv",
                          extras = [species])
    newPipeline.set_head_tasks([headTask])
    newPipeline.set_tail_tasks([tailTask])
    return newPipeline                                   

#---------------------------------------------------------------
# RUN THE PIPELINE
#---------------------------------------------------------------

# make the trimmming pipelines
osjTrimming = makeTrimmingPipeline(pipelineName = "osjTrimming",
                            pathToReads = "data/reads/osj",
                            species = "osj")
osiTrimming = makeTrimmingPipeline(pipelineName = "osiTrimming",
                            pathToReads = "data/reads/osi",
                            species = "osi")
orTrimming = makeTrimmingPipeline(pipelineName = "orTrimming",
                            pathToReads = "data/reads/or",
                            species = "or")
ogTrimming = makeTrimmingPipeline(pipelineName = "ogTrimming",
                            pathToReads = "data/reads/og",
                            species = "og")
obTrimming = makeTrimmingPipeline(pipelineName = "obTrimming",
                            pathToReads = "data/reads/ob",
                            species = "ob")                            

# load the genome into memory
genomeLoad = main_pipeline.transform(task_func = loadGenome_sh,
                                   input = starGenomeGenerate_sh,
                                   filter = suffix("METADATA.csv"),
                                   output = "genomeLoad/METADATA.csv")\
                        .follows(osjTrimming)\
                        .follows(osiTrimming)\
                        .follows(orTrimming)\
                        .follows(ogTrimming)\
                        .follows(obTrimming)\

# run the first step mapping tasks
osjFirstStep = main_pipeline.transform(name = "osjFirstStep",
                                       task_func = firstMapping_sh,
                                       input = osjTrimming,
                                       filter = suffix("cutadapt/METADATA.csv"),
                                       output = "STAR/METADATA.csv",
                                       extras = ["osj"])\
                            .follows(genomeLoad)
osiFirstStep = main_pipeline.transform(name = "osiFirstStep",
                                       task_func = firstMapping_sh,
                                       input = osiTrimming,
                                       filter = suffix("cutadapt/METADATA.csv"),
                                       output = "STAR/METADATA.csv",
                                       extras = ["osi"])\
                            .follows(genomeLoad)
orFirstStep = main_pipeline.transform(name = "orFirstStep",
                                       task_func = firstMapping_sh,
                                       input = orTrimming,
                                       filter = suffix("cutadapt/METADATA.csv"),
                                       output = "STAR/METADATA.csv",
                                       extras = ["or"])\
                            .follows(genomeLoad)
obFirstStep = main_pipeline.transform(name = "obFirstStep",
                                       task_func = firstMapping_sh,
                                       input = obTrimming,
                                       filter = suffix("cutadapt/METADATA.csv"),
                                       output = "STAR/METADATA.csv",
                                       extras = ["ob"])\
                            .follows(genomeLoad)
ogFirstStep = main_pipeline.transform(name = "ogFirstStep",
                                       task_func = firstMapping_sh,
                                       input = ogTrimming,
                                       filter = suffix("cutadapt/METADATA.csv"),
                                       output = "STAR/METADATA.csv",
                                       extras = ["og"])\
                            .follows(genomeLoad)
                             
# unload the genome
genomeUnload = main_pipeline.transform(task_func = unloadGenome_sh,
                                   input = genomeLoad,
                                   filter = suffix("METADATA.csv"),
                                   output = "genomeUnloaded")\
                        .follows(osjFirstStep)\
                        .follows(osiFirstStep)\
                        .follows(orFirstStep)\
                        .follows(ogFirstStep)\
                        .follows(ogFirstStep)\



# run the second step mapping tasks

# print the flowchart
pipeline_printout_graph("ruffus/flowchart." + slurm_jobid + ".pdf", "pdf")

# run the pipeline (disabled for now)
cmdline.run(options, multithread = 8)