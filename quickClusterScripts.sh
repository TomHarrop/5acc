#!/bin/bash

# download genomes
curl https://signon.jgi.doe.gov/signon/create --data-ascii \
	login=thomas.harrop@ird.fr\&password=***REMOVED*** \
	-b cookies -c cookies > /dev/null
	
curl http://genome.jgi.doe.gov/Osativa/download/_JAMO/53112abc49607a1be00559bc/Osativa_204_v7.0.fa.gz \
	-b cookies -c cookies > star/genome/Osativa_204_v7.0.fa.gz

curl http://genome.jgi.doe.gov//Osativa/download/_JAMO/53112ab649607a1be00559b0/Osativa_204_v7.0.gene_exons.gff3.gz \
	-b cookies -c cookies > star/genome/Osativa_204_v7.0.gene_exons.gff3.gz

gunzip star/genome/Osativa_204_v7.0.fa.gz
gunzip star/genome/Osativa_204_v7.0.gene_exons.gff3.gz

rm cookies

# remove Chr9 rRNA 'genes'
sed '/LOC_Os09g01000/d' star/genome/Osativa_204_v7.0.gene_exons.gff3 \
	| sed '/LOC_Os09g00999/d' > star/genome/Osativa_204_v7.0.gene_exons.rRNAremoved.gff3