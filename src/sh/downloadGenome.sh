#!/bin/bash

set -e

# catch email address (-e) and password (-p) for jgi passed to bash script
# if the jgi password has special characters this script won't work

while [ "$1" != "" ]; do
	case $1 in
		-e )	shift
				EMAIL=$1
				;;
		-p )	shift
				PASSWORD=$1
				;;
		* )		echo "Bad input"
				exit 1
	esac
	shift
done

genome_dir='data/genome/os'
genome_url='http://genome.jgi.doe.gov/Osativa/download/_JAMO/5693356c0d87851ee9726b00/Osativa_323_v7.0.fa.gz'
annot_url='http://genome.jgi.doe.gov/Osativa/download/_JAMO/5693356b0d87851ee9726af7/Osativa_323_v7.0.gene_exons.gff3.gz'

genome_file="$(basename $genome_url .fa.gz)"
annotation_file="$(basename $annot_url .gff3.gz)"

# download genomes
if [ ! -d $genome_dir ]; then
	mkdir -p $genome_dir
fi

echo -e "[ "$(date)": Signing on to phytozome at JGI ]"
curl https://signon.jgi.doe.gov/signon/create --data-ascii \
	login="$EMAIL"\&password="$PASSWORD" \
	-b "$genome_dir"/cookies -c "$genome_dir"/cookies > /dev/null

cat <<- _EOF_
	[ $(date): Downloading genome fasta ]
	$genome_url
_EOF_
curl $genome_url -b "$genome_dir"/cookies -c "$genome_dir"/cookies > $genome_dir/$genome_file.fa.gz

# make sure the file was downloaded
if [[ ! -s "$genome_dir"/"$genome_file.fa.gz" ]]; then
	echo -e "[ "$(date)": ERROR: download failed, check password? ]"
	exit 1
fi

cat <<- _EOF_
	[ $(date): Downloading annotation ]
	$annot_url
_EOF_
curl $annot_url -b "$genome_dir"/cookies -c "$genome_dir"/cookies > $genome_dir/$annotation_file.gff3.gz

echo -e "[ "$(date)": Unzipping downloads ]"
gunzip $genome_dir/$genome_file.fa.gz
gunzip $genome_dir/$annotation_file.gff3.gz

rm "$genome_dir"/cookies

# make cuffcomp gtf
echo -e "[ "$(date)": Making GTF file with cuffcompare ]"
cuffcompare -s $genome_dir/$genome_file.fa -CG -r $genome_dir/$annotation_file.gff3 \
	-o $genome_dir/$annotation_file.cuffcomp $genome_dir/$annotation_file.gff3

# remove Chr9 rRNA 'genes'
sed '/LOC_Os09g01000/d' $genome_dir/$annotation_file.cuffcomp.combined.gtf \
	| sed '/LOC_Os09g00999/d' \
	> $genome_dir/gtf_final.tmp

# remove cuffcomp intermediates
rm $genome_dir/*cuffcomp*
mv $genome_dir/gtf_final.tmp $genome_dir/$annotation_file.cuffcomp.rRNAremoved.gtf

# download MSU annotations file
annotationsUrl="ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir/all.locus_brief_info.7.0"
annotationsFile="$(basename $annotationsUrl).tab"
cat <<- _EOF_
	[ $(date): Downloading MSU gene function annotations ]
	$genome_url
_EOF_
curl $annotationsUrl > "$genome_dir"/"$annotationsFile"

# log metadata
cat -t <<- _EOF_ > $genome_dir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cuffcomp version,$(cuffcompare -h 2>&1 | head -n 1)
_EOF_

exit 0