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
 
path_to_dir="harrop@ricebkp.ird.fr:/data/RICE/harrop/"
path_to_tmp="/scratch/tom-$JOB_ID"
path_to_out="/scratch/tom-$JOB_ID-out"

###### Catch arguments

echo -e "\n[ `date`: Preparing to run STAR for library $library_name ] :\nfwd_read_file=$fwd_read_file\nrev_read_file=$rev_read_file\n"
fwd_basename=$(basename $fwd_read_file)
rev_basename=$(basename $rev_read_file)
R1=$path_to_tmp/$fwd_basename
R2=$path_to_tmp/$rev_basename
index_dir=/data/projects/evoreprice/Tom/star/index

###### Create temporary folders on node
mkdir -p $path_to_tmp
mkdir -p $path_to_out/$library_name
echo -e "\n[ `date`: Data transfer master -> node ]\n"
scp -rp $fwd_read_file $R1
scp -rp $rev_read_file $R2

ls -la $path_to_tmp/

###### Run STAR

# ENCODE options
cmd="/usr/local/STAR-2.4.1.c/source/STAR --runThreadN 4 --genomeDir $index_dir --readFilesIn $R1 $R2 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix $path_to_out/$library_name --outSAMtype BAM Unsorted --outSJfilterReads Unique --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outSJfilterIntronMaxVsReadN 5000"

echo -e "\n[ `date`: Running command ]\n$cmd\n"
$cmd

##### Transfer results
echo -e "\n[ `date`: Results transfer node -> master ]\n"
scp -rp $path_to_out $path_to_dir

#### Suppression du repertoire tmp noeud
echo -e "\n[ `date`: Deleting temporary files on node ]\n"
rm -rf $path_to_tmp
rm -rf $path_to_out
