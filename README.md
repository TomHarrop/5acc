## five_accessions

Analysis for https://github.com/TomHarrop/ird-5acc-paper

### Requirements

[`singularity`](https://singularity.lbl.gov)  
[`snakemake`](https://snakemake.readthedocs.io) ≥ 4.7.0

### Reproduce the analysis

`snakemake --use-singularity --cores={threads} --resources mem_gb={ram_limit}`

### Software environment

The software used in the analysis is listed below. A [`singularity`](singularity.lbl.gov) container for this pipeline is hosted at [shub://TomHarrop/singularity-containers:five-accessions](https://www.singularity-hub.org/collections/996).

[`snakemake`](https://snakemake.readthedocs.io) will pull the container using the `--use-singularity` flag, and all analysis will run with software installed in the container.

- `bbmap` 38.00
- `bedtools` 2.26.0
- `cuffcompare` 2.2.1
- `STAR` 2.5.4b
- `wgsim` 0.3.1-r13
- `python` 3.6.5, with packages:
    + `cutadapt` 1.16
    + `HTSeq` 0.9.1
- `R` 3.4.4, with packages
    + `data.table` 1.11.2
    + `DESeq2` 1.18.1
    + `GenomicRanges` 1.30.3
    + `rtracklayer`  1.30.3
    + `valr` 0.4.0  

### Input data files

The following files are not distributed with the workflow, and must be in a `data` directory under the current working directory.

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

