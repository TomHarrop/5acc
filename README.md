## five_accessions

Analysis for https://github.com/TomHarrop/ird-5acc-paper

### Requirements

`bbmap` script `repair.sh`  
`bedtools` ≥ 2.26.0   
`cuffcompare` ≥ 2.2.1  
`STAR` ≥ 2.5.3a  
`wgsim` ≥ 0.3.2  

#### `python` packages

`cutadapt` ≥ 1.16  
`snakemake` ≥ 4.7.0  

#### `R` packages

`data.table` ≥ 1.10.4-3  
`GenomicRanges` ≥ 1.30.3  
`rtracklayer` ≥ 1.38.3  
`valr` ≥ 0.4.0  

### Input data files

`data/genome/os/Osativa_323_v7.0.fa`  
`data/genome/os/Osativa_323_v7.0.gene_exons.gff3`  
`data/reads/{sample}_R1.fastq.gz`  
`data/reads/{sample}_R2.fastq.gz`  

### Pipeline

![](dag/dag.svg)

### Run the analysis

`snakemake --cores={threads} --resources mem_gb={ram_limit}`
