## five_accessions

Analysis for https://github.com/TomHarrop/ird-5acc-paper

### Requirements

`bbmap` script `repair.sh`  
`cuffcompare` ≥ 2.2.1  
`cutadapt` ≥ 1.16
`snakemake` ≥ 4.7.0  
`STAR` ≥ 2.5.3a  

### Input data files

`data/genome/os/Osativa_323_v7.0.fa`  
`data/genome/os/Osativa_323_v7.0.gene_exons.gff3`  
`data/reads/{sample}_R1.fastq.gz`  
`data/reads/{sample}_R2.fastq.gz`  

### Pipeline

![](text/dag.svg)

### Run the analysis

`snakemake --cores={threads} --resources mem_gb={ram_limit}`