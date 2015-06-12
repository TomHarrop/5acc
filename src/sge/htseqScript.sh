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
#$ -N htseq

############################################################
 
###### Catch arguments

GTF="/data/projects/evoreprice/Tom/star/genome/Osativa_204_v7.0.gene_exons.rRNAremoved.gff3"
echo -e "[ `date`: Preparing to run htseq-count for library $library_name ]\nParamaters passed from wrapper:\nlibrary_name:\t\t$library_name\nbamfile:\t\t$bamfile"
echo -e "Using annotation:\t$GTF"

bamfile_base="$(basename $bamfile)"
GTF_base="$(basename $GTF)"

###### trap interruptions (not working, SGE sends SIGKILL???)

clean_up() {
	# remove temp files before exit
	echo -e "[ `date`: script aborted! Deleting temporary files ]"
	rm -rf $path_to_tmp
	rm -rf $path_to_out
	exit
}

trap clean_up SIGHUP SIGINT SIGTERM

###### Create temporary folders on node

path_to_dir="/data/projects/evoreprice/Tom/htseq-count"
path_to_tmp="/scratch/tom.$library_name.htseq.$JOB_ID"
path_to_out="/scratch/tom.$library_name.htseq.$JOB_ID.out"

mkdir -p $path_to_tmp
mkdir -p $path_to_out

GTF_tmp="$path_to_tmp/$GTF_base"
bamfile_tmp="$path_to_tmp/$bamfile_base"

echo -e "[ `date`: Data transfer master -> node ]"
scp $bamfile $bamfile_tmp
scp $GTF $GTF_tmp

echo -e "[ `date`: Finished copying ]"
ls -lh /
ls -lh $path_to_tmp
find $path_to_tmp/*

###### Run program

cmd="python -m HTSeq.scripts.count -f bam -s reverse -i Parent $bamfile_tmp $GTF_tmp > $path_to_out/$library_name.htseq-count"
echo -e "[ `date`: Running command ]\n$cmd"

python -m HTSeq.scripts.count -f bam -s reverse -i Parent $bamfile_tmp $GTF_tmp > $path_to_out/$library_name.htseq-count

##### Transfer results
echo -e "[ `date`: Results transfer node -> master ]"
scp -rp $path_to_out $path_to_dir

#### Suppression du repertoire tmp noeud
echo -e "[ `date`: Tranfer finished. Deleting temporary files on node ]"
rm -rf $path_to_tmp
rm -rf $path_to_out

echo -e "[ `date`: Done ]\n"
