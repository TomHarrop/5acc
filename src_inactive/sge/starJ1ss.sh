#!/bin/bash

############      SGE CONFIGURATION      ###################
# join error 
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
#$ -q highmem.q

# Nom du job
#$ -N star

# Resources
#$ -pe smp 4

############################################################
 
export LANG="en_AU.UTF-8"
 
###### Run STAR

gff="/data/projects/evoreprice/Tom/star/genome/Osativa_204_v7.0.gene_exons.gff3"
fasta="/data/projects/evoreprice/Tom/star/genome/Osativa_204_v7.0.fa"
index_dir="/data/projects/evoreprice/Tom/star/index/"
R1="/data/projects/evoreprice/Tom/subsample/J1/R1.fastq"
R2="/data/projects/evoreprice/Tom/subsample/J1/R2.fastq"

# default options
echo -e "[ $(date): Running STAR ]"
/home/harrop/dl/STAR/source/STAR --runThreadN 4 --genomeDir "$index_dir" \
	--readFilesIn "$R1" "$R2" --outFileNamePrefix /data/projects/evoreprice/Tom/subsample/J1/STAR/J1ss. \
	--outSAMtype BAM Unsorted
echo -e "[ $(date): Done ]"
exit 0
