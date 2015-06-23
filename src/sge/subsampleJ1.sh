#!/bin/bash

############      SGE CONFIGURATION      ###################
# Join error
#$ -j y

# Which shell
#$ -S /bin/bash 

# Email pour suivre l'execution 
#$ -M thomas.harrop@ird.fr

# Type de massage que l'on reçoit par mail
#    -  (b) un message au demarrage
#    -  (e) a la fin
#    -  (a)  en cas d'abandon
#$ -m bea 

# Queue que l'on veut utiliser
#$ -q bioinfo.q

# Nom du job
#$ -N subsamp

############################################################

export LANG="en_AU.UTF-8"

path_to_tmp="/scratch/tom-subsample-"$JOB_ID""

mkdir -p "$path_to_tmp"

echo -e "[ $(date): gunzip fastq files ]"

gunzip -c \
	/data/projects/evoreprice/Data/OJ/Trimmed/OJ-cutadapt_cluster-20150522/J1_AGTCAA_L002_R1.fastq.fastq.gz \
	> "$path_to_tmp"/R1.fastq
gunzip -c \
	/data/projects/evoreprice/Data/OJ/Trimmed/OJ-cutadapt_cluster-20150522/J1_AGTCAA_L002_R2.fastq.fastq.gz \
	> "$path_to_tmp"/R2.fastq

echo -e "[ $(date): done ]"
ls -lh "$path_to_tmp"

echo -e "[ $(date): subsampling ]"
/usr/local/bin/python /home/harrop/scripts/subsampler_PE_anders.py \
	0.002449463 "$path_to_tmp"/R1.fastq "$path_to_tmp"/R2.fastq \
	/data/projects/evoreprice/Tom/subsample/J1/R1.fastq \
	/data/projects/evoreprice/Tom/subsample/J1/R2.fastq

echo -e "[ $(date): done, cleaning up ]"
rm -r "$path_to_tmp"
exit 0
