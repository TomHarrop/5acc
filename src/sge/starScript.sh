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
 
path_to_dir="/data/projects/evoreprice/Tom"
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
echo -e "\n[ `date`: Data transfer master -> node ]"
scp -rp $fwd_read_file $R1
scp -rp $rev_read_file $R2

ls -la $path_to_tmp/

###### Run STAR

gff=$path_to_tmp/genome/Osativa_204_v7.0.gene_exons.rRNAremoved.gff3
fasta=$path_to_tmp/genome/Osativa_204_v7.0.fa

# default options
echo -e "\n[ `date`: Running STAR ]"
/home/harrop/dl/STAR/source/STAR --runThreadN 4 --genomeDir $index_dir \
	--readFilesIn $R1 $R2 --outFileNamePrefix $path_to_out/$library_name. \
	--outSAMtype BAM Unsorted

##### Transfer results
echo -e "\n[ `date`: Results transfer node -> master ]"
rcp -rp $path_to_out $path_to_dir

#### Suppression du repertoire tmp noeud
echo -e "\n[ `date`: Deleting temporary files on node ]"
rm -rf $path_to_tmp
rm -rf $path_to_out
