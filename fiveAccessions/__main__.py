#!/usr/bin/python3
# -*- coding: utf-8 -*-
#

#########################
# 5 accessions pipeline #
#########################

import functions
import ruffus
import os
import datetime

###################
# SETUP FUNCTIONS #
###################


# download rice genome
def downloadGenome_sh(outputFiles, jgi_logon, jgi_password):
    jobScript = 'src/sh/downloadGenome.sh'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'downloadGenome_sh'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=['-e', jgi_logon, '-p', jgi_password])
    functions.print_job_submission(job_name, job_id)


# generate STAR index for OS
def starGenomeGenerate_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/starGenomeGenerate.sh'
    ntasks = '1'
    cpus_per_task = '8'
    job_name = 'stargg'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)


# calculate feature lengths from GTF
def calculate_feature_lengths_R(inputFiles, outputFiles):
    jobScript = 'src/R/calculate_feature_lengths.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = 'feature_lengths'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)


# shuffle the GTF
def shuffle(inputFiles, outputFiles):
    jobScript = 'src/sh/shuffle'
    ntasks = '1'
    cpus_per_task = '8'
    job_name = 'shuffle'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)

##################################
# MAPPING AND TRIMMING FUNCTIONS #
##################################


# define reads
def defineReads(outputFiles, species):
    pathToReads = "data/reads/" + species
    assert os.path.isdir(pathToReads), ("1 " + outputFiles + " 2 " +
                                        pathToReads + " 3 " + species +
                                        "Error: reads folder " + pathToReads +
                                        " missing")
    readFiles = os.listdir(pathToReads)
    print("[", datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
          ": Using reads in folder " + pathToReads + " ]")
    for fileName in readFiles:
        qualName = pathToReads + '/' + fileName
        assert os.path.isfile(qualName), ("Error: read file " + qualName +
                                          " missing")
        print(qualName)
    functions.touch(outputFiles)


# trim reads
def cutadapt_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/cutadapt.sh'
    ntasks = '2'
    cpus_per_task = '1'
    job_name = species + '_cutadapt'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=[species])
    functions.print_job_submission(job_name, job_id)


# load genome into memory
def loadGenome_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/loadGenome.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'loadGenome'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)


# first mapping step (shared memory across jobs)
def firstMapping_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/firstMapping.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = species + '_firstMapping'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=[species])
    functions.print_job_submission(job_name, job_id)


# second mapping step (can't use shared memory with 2-pass mapping)
def secondMapping_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/secondMapping.sh'
    ntasks = '1'
    cpus_per_task = '8'
    job_name = species + '_secondMapping'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=[species])
    functions.print_job_submission(job_name, job_id)


# unload genome from memory
def unloadGenome_sh(inputFiles, outputFiles):
    jobScript = 'src/sh/unloadGenome.sh'
    ntasks = '1'
    cpus_per_task = '7'
    job_name = 'unloadGenome'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)


# sort bamfiles


# compress unmapped read files
def compressUnmappedReads_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/compressUnmappedReads.sh'
    ntasks = '1'
    cpus_per_task = '8'
    job_name = species + '_gzip'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=[species])
    functions.print_job_submission(job_name, job_id)


# remap reads
def remap_sh(inputFiles, outputFiles, species):
    jobScript = 'src/sh/remap.sh'
    ntasks = '1'
    cpus_per_task = '8'
    job_name = species + '_remap'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=[species])
    functions.print_job_submission(job_name, job_id)


# count reads in shuffled GFF
def htseq_shuffle(inputFiles, outputFiles, species):
    jobScript = 'src/sh/htseq_shuffle'
    ntasks = '4'
    cpus_per_task = '1'
    job_name = species + '_count'
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name,
                                  extras=["-s", species])
    functions.print_job_submission(job_name, job_id)


#################################
# DOWNSTREAM ANALYSIS FUNCTIONS #
#################################


# parse mapping stats
def parseStarStats_R(inputFiles, outputFile):
    jobScript = 'src/R/parseStarStats.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = "parseStarStats"
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)


# run DESeq2
def deseq2_R(inputFiles, outputFiles):
    jobScript = 'src/R/DESeq2.R'
    ntasks = '1'
    cpus_per_task = '1'
    job_name = "deseq2_R"
    job_id = functions.submit_job(jobScript, ntasks, cpus_per_task, job_name)
    functions.print_job_submission(job_name, job_id)


# calculate TPM
def tpm_R(inputFiles, outputFiles):
    pass


# perform QC checks on DESeq2 output
def deseqQC_R(inputFiles, outputFiles):
    pass


# calculate expression cutoffs
def cutoffs_R(inputFiles, outputFiles, species):
    pass


##############################
# IMPLEMENT POST-RUN BACKUP? #
##############################

##########################
# CONSTRUCT THE PIPELINE #
##########################

def main():
    # get default pipeline
    main_pipeline = ruffus.Pipeline.pipelines["main"]

    # catch jgi logon and password from cli
    parser = ruffus.cmdline.get_argparse(description='5 accessions analysis '
                                                     'pipeline.')
    parser.add_argument('--email', '-e',
                        help='Logon email address for JGI',
                        type=str,
                        dest='jgi_logon')
    parser.add_argument('--password', '-p',
                        help='JGI password',
                        type=str,
                        dest='jgi_password')
    options = parser.parse_args()
    jgi_logon = options.jgi_logon
    jgi_password = options.jgi_password

    # download the genome
    os_genome = main_pipeline.originate(downloadGenome_sh,
                                        "data/genome/os/METADATA.csv",
                                        jgi_logon, jgi_password)

    # create the STAR index
    os_index = main_pipeline.transform(
        task_func=starGenomeGenerate_sh,
        input=os_genome,
        filter=ruffus.suffix("data/genome/os/METADATA.csv"),
        output="output/star-index/METADATA.csv")

    # calculate feature lengths
    feature_lengths = main_pipeline.transform(
        task_func=calculate_feature_lengths_R,
        input=os_genome,
        filter=ruffus.suffix("data/genome/os/METADATA.csv"),
        output="output/tpm/SessionInfo.calculate_feature_lengths.txt")

    # shuffle the GTF
    shuffled_gff = main_pipeline.transform(
        task_func=shuffle,
        input=feature_lengths,
        filter=ruffus.suffix("tpm/SessionInfo.calculate_feature_lengths.txt"),
        output="shuffle/METADATA.csv")

    # define the reads
    species = ["osj", "osi", "or", "ob", "og"]
    for spec in species:
        main_pipeline.originate(name=spec + "Reads",
                                task_func=defineReads,
                                output="ruffus/" + spec + ".reads",
                                extras=[spec])

    # run cutadapt
    trimming = main_pipeline.transform(
        task_func=cutadapt_sh,
        input=ruffus.output_from(list(n + "Reads" for n in species)),
        filter=ruffus.regex(r"ruffus/(.*).reads"),
        output=r"output/\1/cutadapt/METADATA.csv",
        extras=[r"\1"])

    # load the genome into memory
    genomeLoad = main_pipeline.transform(task_func=loadGenome_sh,
                                         input=starGenomeGenerate_sh,
                                         filter=ruffus.suffix("METADATA.csv"),
                                         output="genomeLoad/METADATA.csv")\
                              .follows(trimming)

    # run the first step mapping tasks
    firstStep = main_pipeline.transform(
        task_func=firstMapping_sh,
        input=trimming,
        filter=ruffus.regex(r"output/(.*)/cutadapt/METADATA.csv"),
        output=r"output/\1/STAR/step1/METADATA.csv",
        extras=[r"\1"])\
        .follows(genomeLoad)

    # unload the genome
    genomeUnload = main_pipeline.transform(
        task_func=unloadGenome_sh,
        input=genomeLoad,
        filter=ruffus.suffix("METADATA.csv"),
        output="genomeUnloaded")\
        .follows(firstStep)

    # run the second step mapping tasks
    secondStep = main_pipeline.transform(
        task_func=secondMapping_sh,
        input=firstStep,
        filter=ruffus.regex(r"output/(.*)/STAR/step1/METADATA.csv"),
        output=r"output/\1/STAR/METADATA.csv",
        extras=[r"\1"])\
        .follows(genomeUnload)

    # parse stats from the mapping
    parsed_stats = main_pipeline.merge(
        task_func=parseStarStats_R,
        input=secondStep,
        output="output/mappingStats/starLogs.Rds")

    # count reads in shuffled GFF
    shuffled_reads = main_pipeline.transform(
        task_func=htseq_shuffle,
        input=secondStep,
        add_inputs=ruffus.add_inputs(shuffled_gff),
        filter=ruffus.formatter(),
        output=["output/shuffle/htseq/{subdir[0][1]}/METADATA.csv"],
        extras=["{subdir[0][1]}"])

    # run deseq2 for basic comparisons and QC tasks
    deseq2 = main_pipeline.merge(task_func=deseq2_R,
                                 input=secondStep,
                                 output="output/deseq2/SessionInfo.txt")

    # convert transformed read counts to TPM
    tpm = main_pipeline.merge(
        task_func=tpm_R,
        input=[parsed_stats, deseq2, feature_lengths],
        output=["output/tpm/SessionInfo.tpm.txt"])

    # compress unmapped read files
    compress = main_pipeline.transform(
        task_func=compressUnmappedReads_sh,
        input=secondStep,
        filter=ruffus.regex(r"output/(.*)/STAR/METADATA.csv"),
        output=r"output/\1/STAR/compressStats.txt",
        extras=[r"\1"])

    # remap unmapped read files
    remap = main_pipeline.transform(
        task_func=remap_sh,
        input=compress,
        filter=ruffus.regex(r"output/(.*)/STAR/compressStats.txt"),
        output=r"output/\1/STAR/remap/METADATA.csv",
        extras=[r"\1"])

    # run QC on deseq2 output
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

    # COMMAND-LINE OPTIONS

    # print the flowchart
    ruffus.pipeline_printout_graph("ruffus/flowchart.pdf", "pdf",
                                   pipeline_name="5 accessions analysis "
                                                 "pipeline")

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

if __name__ == "__main__":
    main()
