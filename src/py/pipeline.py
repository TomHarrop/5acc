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

#---------------------------------------------------------------
# MAPPING AND TRIMMING TASKS
#---------------------------------------------------------------

#---------------------------------------------------------------
# define reads
#
def defineReads(outputFiles, species):
    pathToReads = "data/reads/" + species
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
    job_name = species + '_cutadapt'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = species)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")
      
#---------------------------------------------------------------
# load genome into memory
#
def loadGenome_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/loadGenome.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'loadGenome'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# first mapping step (shared memory across jobs)
#
def firstMapping_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/firstMapping.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = species + '_firstMapping'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = species)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# second mapping step (can't use shared memory with 2-pass mapping)
#
def secondMapping_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/secondMapping.sh'
    ntasks = '1'
    cpus_per_task = '8'
    job_name = species + '_secondMapping'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = species)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# unload genome from memory
#
def unloadGenome_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/unloadGenome.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'unloadGenome'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# sort bamfiles
#


#---------------------------------------------------------------
# compress unmapped read files
#
def compressUnmappedReads_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/compressUnmappedReads.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = species + '_gzip'
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name, extras = species)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# remap reads
#

#---------------------------------------------------------------
# DOWNSTREAM ANALYSIS
#---------------------------------------------------------------

#---------------------------------------------------------------
# parse mapping stats
#
def parseStarStats_R(inputFiles, outputFile):
    jobScript= 'src/R/parseStarStats.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = "parseStarStats"
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# run DESeq2
#
def deseq2_R(inputFiles, outputFiles):
    jobScript= 'src/R/DESeq2.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = "deseq2_R"
    jobId = submit_job(jobScript, ntasks, cpus_per_task, job_name)
    print("[", print_now(), ": Job " + job_name + " run with JobID " + jobId + " ]")

#---------------------------------------------------------------
# perform QC checks on DESeq2 output
#
def deseqQC_R(inputFiles, outputFiles):
    pass

#---------------------------------------------------------------
# calculate expression cutoffs
#
def cutoffs_R(inputFiles, outputFiles, species):
    pass

    
#---------------------------------------------------------------
# RUN THE PIPELINE
#---------------------------------------------------------------

# get default pipeline
main_pipeline = Pipeline.pipelines["main"]

# download the genome
os_genome = main_pipeline.originate(downloadGenome_sh, "data/genome/os/METADATA.csv", jgiLogon, jgiPassword)

# create the STAR index
os_index = main_pipeline.transform(task_func = starGenomeGenerate_sh,
                                   input = downloadGenome_sh,
                                   filter = suffix("data/genome/os/METADATA.csv"),
                                   output = "output/star-index/METADATA.csv")

# define the reads
species = ["osj", "osi", "or", "ob", "og"]
for spec in species:
    main_pipeline.originate(name = spec + "Reads",
                            task_func = defineReads,
                            output = "ruffus/" + spec + ".reads",
                            extras = [spec])

# run cutadapt
trimming = main_pipeline.transform(task_func = cutadapt_sh,
                                   input = output_from(list(n + "Reads" for n in species)),
                                   filter = regex(r"ruffus/(.*).reads"),
                                    output = r"output/\1/cutadapt/METADATA.csv",
                                    extras = [r"\1"])

# load the genome into memory
genomeLoad = main_pipeline.transform(task_func = loadGenome_sh,
                                   input = starGenomeGenerate_sh,
                                   filter = suffix("METADATA.csv"),
                                   output = "genomeLoad/METADATA.csv")\
                        .follows(trimming)

# run the first step mapping tasks
firstStep = main_pipeline.transform(task_func = firstMapping_sh,
                                       input = trimming,
                                       filter = regex(r"output/(.*)/cutadapt/METADATA.csv"),
                                       output = r"output/\1/STAR/step1/METADATA.csv",
                                       extras = [r"\1"])\
                            .follows(genomeLoad)
                            
# unload the genome                            
genomeUnload = main_pipeline.transform(task_func = unloadGenome_sh,
                                   input = genomeLoad,
                                   filter = suffix("METADATA.csv"),
                                   output = "genomeUnloaded")\
                        .follows(firstStep)
                        
# run the second step mapping tasks    
secondStep = main_pipeline.transform(task_func = secondMapping_sh,
                                       input = firstStep,
                                       filter = regex(r"output/(.*)/STAR/step1/METADATA.csv"),
                                       output = r"output/\1/STAR/METADATA.csv",
                                       extras = [r"\1"])\
                            .follows(genomeUnload)

# parse stats from the mapping
parseStats = main_pipeline.merge(task_func = parseStarStats_R,
                                 input = secondStep,
                                 output = "output/mappingStats/starLogs.Rds")

# run deseq2 for basic comparisons and QC tasks
deseq2 = main_pipeline.merge(task_func = deseq2_R,
                                 input = secondStep,
                                 output = "output/deseq2/SessionInfo.txt")

# compress unmapped read files
compress = main_pipeline.transform(task_func = compressUnmappedReads_sh,
                                   input = secondStep,
                                   filter = regex(r"output/(.*)/STAR/METADATA.csv"),
                                    output = r"output/\1/STAR/compressStats.txt",
                                    extras = [r"\1"])

## run QC on deseq2 output
#deseqQC = main_pipeline.transform(task_func = deseqQC_R,
#                                  input = deseq2,
#                                  filter = suffix("SessionInfo.txt"),
#                                  output = "someFileHere.csv")

# calculate cutoffs
#deseq2 = main_pipeline.transform(task_func = cutoffs_R,
#                                 input = secondStep,
#                                 filter = regex(r"output/(.*)/STAR/METADATA.csv"),
#                                 output = r"output/\1/cutoffs/SessionInfo.txt",
#                                 extras = [r"\1"])

#---------------------------------------------------------------
# COMMAND-LINE OPTIONS
#---------------------------------------------------------------

# print the flowchart
pipeline_printout_graph("ruffus/flowchart." + slurm_jobid + ".pdf", "pdf",
                        pipeline_name = "5accessions analysis pipeline")

# run the pipeline (disabled for now)
cmdline.run(options, multithread = 8)