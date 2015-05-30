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
#$ -q bioinfo.q

# Nom du job
#$ -N starGen

# Resources
#$ -pe smp 4

############################################################
 
path_to_dir="/data/projects/evoreprice/Tom"
path_to_tmp="/scratch/tom-$JOB_ID"
path_to_out="/scratch/tom-$JOB_ID-out"
 
###### Create temporary folders on node
mkdir -p $path_to_tmp
mkdir -p $path_to_out/star/index
echo -e "\n[ `date`: Data transfer master -> node ]\n"
scp -rp $path_to_dir/star/genome $path_to_tmp/
ls -la $path_to_tmp/
ls -la $path_to_tmp/genome/

###### Run STAR

gff=$path_to_tmp/genome/Osativa_204_v7.0.gene_exons.rRNAremoved.gff3
fasta=$path_to_tmp/genome/Osativa_204_v7.0.fa

cmd="/usr/local/STAR-2.4.1.c/source/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $path_to_out/star/index --genomeFastaFiles $fasta --sjdbGTFfile $gff --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 109"
echo -e "\n[ `date`: Running command ]\n$cmd\n"
$cmd

##### Transfer results
echo -e "\n[ `date`: Results transfer node -> master ]\n"
rcp -rp $path_to_out $path_to_dir

#### Suppression du repertoire tmp noeud
echo -e "\n[ `date`: Deleting temporary files on node ]\n"
rm -rf $path_to_tmp
rm -rf $path_to_out
