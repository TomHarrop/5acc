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
`HTSeq` ≥ 0.9.1  
`snakemake` ≥ 4.7.0  

#### `R` packages

`data.table` ≥ 1.10.4-3  
`DESeq2` ≥ 1.18.1  
`GenomicRanges` ≥ 1.30.3  
`rtracklayer` ≥ 1.38.3  
`valr` ≥ 0.4.0  

### Input data files

- Raw reads:
    + `data/reads/{sample}_R1.fastq.gz`
    + `data/reads/{sample}_R2.fastq.gz`
- From Phytozome:
    + `data/genome/os/Osativa_323_v7.0.fa`
    + `data/genome/os/Osativa_323_v7.0.gene_exons.gff3`
- From http://rapdb.dna.affrc.go.jp/download/archive:
    + `data/genome/os/irgsp1_rRNA_tRNA.gff`
- From ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir:
    + `data/genome/os/rice_osa1r7_rm.gff3`
- From ftp://ftp.plantbiology.msu.edu/pub/data/TIGR_Plant_Repeats/TIGR_Oryza_Repeats.v3.3:
    + `data/genome/os/TIGR_Oryza_Repeats.v3.3_0_0.fsa`
- From ftp://ftp.plantbiology.msu.edu/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_7.0/all.dir:
    + `data/genome/os/osa.gff3`

### Pipeline

![](dag/dag.svg)

### Run the analysis

`snakemake --cores={threads} --resources mem_gb={ram_limit}`
