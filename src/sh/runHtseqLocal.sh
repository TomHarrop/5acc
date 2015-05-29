#!/bin/bash

# catch paramaters

bamfile=$1
library_name=$2
echo -e "[ `date`: Preparing to run htseq-count for library $library_name ]\nParamaters passed from wrapper:\nlibrary_name:\t\t$library_name\nbamfile:\t\t$bamfile"

bamfile_base=$(basename $bamfile)

# make local directories

local_dir="/home/tom/Desktop/5accessions/output/htseq`date +%F`"

if [ ! -d $local_dir ]; then
	echo "[`date`: making local_dir $local_dir ]"
	mkdir -p $local_dir
else
	echo "[`date`: local_dir $local_dir already exists... skipping ]"
fi

# copy files to local storage

bamfile_local="$local_dir/$bamfile_base"


if [! -e $bamfile_local ]; then
	echo -e "[ `date`: copying bamfile to local_dir ]"
	scp ***REMOVED***:$bamfile $bamfile_local
	echo -e "[ `date`: Finished copying ]"
fi

# run htseq-count

GTF="/home/tom/Desktop/5accessions/data/genome/Osativa_204_v7.0.gene_exons.rRNAremoved.gff3"
error="$local_dir/$library_name.error"
results="$local_dir/$library_name.htseq"
cmd="srun -c 1 -J ht$library_name -e $error -o $results htseq-count -f bam -s reverse -i Parent $bamfile_local $GTF &"

echo -e "[ `date`: Preparing to run htseq-count for library $library_name:]\nbamfile:\t$bamfile_local\nGTF:\t\t$GTF\ncmd:\t\t$cmd"

srun -c 1 -J ht$library_name -e $error -o $results htseq-count -f bam -s reverse -i Parent $bamfile_local $GTF &

wait

# output results to bioinfo-master

echo -e "[ `date`: Transferring results to bioinfo-master ]"

scp -p $error ***REMOVED***:/data/projects/evoreprice/Tom/htseq-count
scp -p $results ***REMOVED***:/data/projects/evoreprice/Tom/htseq-count

echo -e "[ `date`: Done ]"