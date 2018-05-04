#!/usr/bin/env python3

import os


#############
# FUNCTIONS #
#############

def FindAllFastqFiles(read_dir):
    '''Return a list of paths to fastq_files in read_dir'''
    read_dir_files = [(dirpath, filenames) for (dirpath, dirnames, filenames)
                      in os.walk(read_dir, followlinks=True)]
    fastq_files = []
    for dirpath, filenames in read_dir_files:
        for filename in filenames:
            if 'fastq.gz' in filename:
                fastq_files.append(os.path.join(dirpath, filename))
    return(fastq_files)


def FindInputReads(wildcards):
    # globals
    stage_rep_to_sample_number = {
        'PBM_1': '1',
        'PBM_2': '2',
        'PBM_3': '3',
        'SM_1': '4',
        'SM_2': '5',
        'SM_3': '6'}
    species_code_to_letter = {
        'osj': 'J',
        'osi': 'I',
        'ob': 'B',
        'or': 'R',
        'og': 'G'}
    my_stage = wildcards.stage
    my_species = wildcards.species
    my_rep = wildcards.rep
    my_sample_number = stage_rep_to_sample_number[
        '_'.join([my_stage, my_rep])]
    my_letter = species_code_to_letter[my_species]
    my_fastq_name = ''.join([my_letter, my_sample_number])
    my_r1 = [x for x in all_fastq_files
             if('R1' in os.path.basename(x)
                and os.path.basename(x).startswith(my_fastq_name))][0]
    my_r2 = [x for x in all_fastq_files
             if('R2' in os.path.basename(x)
                and os.path.basename(x).startswith(my_fastq_name))][0]
    return({'r1': my_r1, 'r2': my_r2})


###########
# GLOBALS #
###########

os_genome = 'data/genome/os/Osativa_323_v7.0.fa'
os_gff = 'data/genome/os/Osativa_323_v7.0.gene_exons.gff3'
os_gtf = 'output/010_data/Osativa_323_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf'
read_dir = 'data/reads'
all_species = ['osj',
               'osi',
               'ob',
               'or',
               'og']
all_reps = ['1', '2', '3']
all_stages = ['PBM', 'SM']

#########
# SETUP #
#########

all_fastq_files = FindAllFastqFiles(read_dir)

#########
# RULES #
#########

rule target:
    input:
        expand(('output/030_mapping/star-pass2/{species}/'
                '{stage}_{rep}.ReadsPerGene.out.tab'),
               species=all_species,
               stage=all_stages,
               rep=all_reps)

# 030 map
rule second_mapping:
    input:
        r1 = 'output/021_trim-reads/{species}/{stage}_{rep}.r1.fastq.gz',
        r2 = 'output/021_trim-reads/{species}/{stage}_{rep}.r2.fastq.gz',
        genome = 'output/010_data/star-index/SA',
        sjs = expand(('output/030_mapping/star-pass1/{{species}}/'
                      '{stage}_{rep}.SJ.out.tab'),
                     stage=all_stages,
                     rep=all_reps)
    output:
        ('output/030_mapping/star-pass2/{species}/'
         '{stage}_{rep}.ReadsPerGene.out.tab'),
        ('output/030_mapping/star-pass2/{species}/'
         '{stage}_{rep}.Aligned.out.bam')
    threads:
        20
    resources:
        mem_gb = 40
    params:
        genome_dir = 'output/010_data/star-index',
        prefix = 'output/030_mapping/star-pass2/{species}/{stage}_{rep}.'
    log:
        'output/000_logs/030_mapping/{species}-{stage}-{rep}_pass2.log'
    shell:
        'STAR '
        '--sjdbFileChrStartEnd {input.sjs} '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--readFilesCommand zcat '
        '--outSAMtype BAM Unsorted '
        '--quantMode GeneCounts '
        '--outBAMcompression 10 '
        '--outReadsUnmapped Fastx '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'


rule first_mapping:
    input:
        r1 = 'output/021_trim-reads/{species}/{stage}_{rep}.r1.fastq.gz',
        r2 = 'output/021_trim-reads/{species}/{stage}_{rep}.r2.fastq.gz',
        genome = 'output/010_data/star-index/SA'
    output:
        'output/030_mapping/star-pass1/{species}/{stage}_{rep}.SJ.out.tab'
    threads:
        50
    resources:
        mem_gb = 40
    params:
        genome_dir = 'output/010_data/star-index',
        prefix = 'output/030_mapping/star-pass1/{species}/{stage}_{rep}.'
    log:
        'output/000_logs/030_mapping/{species}-{stage}-{rep}_pass1.log'
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--readFilesCommand zcat '
        '--outSAMtype None '
        '--readFilesIn {input.r1} {input.r2} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'


# 020 trim reads
rule gzip:
    input:
        'output/021_trim-reads/{species}/{stage}_{rep}.r{n}.fastq'
    output:
        'output/021_trim-reads/{species}/{stage}_{rep}.r{n}.fastq.gz'
    threads:
        1
    shell:
        'cat {input} | gzip -9 > {output}'

rule cutadapt:
    input:
        r1 = 'output/020_repair/{species}/{stage}_{rep}.r1.fastq',
        r2 = 'output/020_repair/{species}/{stage}_{rep}.r2.fastq'
    output:
        r1 = temp('output/021_trim-reads/{species}/{stage}_{rep}.r1.fastq'),
        r2 = temp('output/021_trim-reads/{species}/{stage}_{rep}.r2.fastq')
    threads:
        1
    log:
        'output/000_logs/021_trim-reads/{species}_{stage}_{rep}.log'
    shell:
        'cutadapt '
        '-a \'TruSeq_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTC\' '
        '-A \'Illumina_single_end=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG\' '
        '--quality-cutoff=20 '
        '--minimum-length=50 '
        '--output={output.r1} '
        '--paired-output={output.r2} '
        '{input.r1} {input.r2} '
        '&> {log}'


rule repair:
    input:
        unpack(FindInputReads)
    output:
        r1 = temp('output/020_repair/{species}/{stage}_{rep}.r1.fastq'),
        r2 = temp('output/020_repair/{species}/{stage}_{rep}.r2.fastq')
    threads:
        1
    resources:
        mem_gb = 50
    log:
        'output/000_logs/020_repair/{species}_{stage}_{rep}.log'
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'zl=9 '
        'repair=t '
        '-Xmx{resources.mem_gb}g '
        '2> {log}'


# 010 prepare data
# rule shuffle_gtf:
#     input:
#         os_gff_file = "data/genome/os/Osativa_323_v7.0.gene_exons.gff3",
#         os_gtf_file = "output/010_data/Osativa_323_v7.0.gene_exons.cuffcomp.rRNAremoved.gtf",
#         seqlengths_file = "output/010_data/star-index/chrNameLength.txt",
#         irgsp_gff_file = "data/genome/os/irgsp1_rRNA_tRNA.gff",
#         osa1r7_gff_file = "data/genome/os/rice_osa1r7_rm.gff3",
#         osa1_mirbase_gff_file = "data/genome/os/osa.gff3",
#         tigr_repeats_fa = "data/genome/os/TIGR_Oryza_Repeats.v3.3_0_0.fsa"
#     params:
#         star_index_dir = "output/010_data/star-index"


rule generate_genome:
    input:
        os_genome = os_genome,
        os_gtf = os_gtf
    output:
        'output/010_data/star-index/SA'
    params:
        outdir = 'output/010_data/star-index'
    threads:
        50
    log:
        'output/000_logs/010_prepare-data/generate_genome.log'
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input.os_genome} '
        '--sjdbGTFfile {input.os_gtf} '
        '--sjdbGTFtagExonParentTranscript oId '
        '--sjdbGTFtagExonParentGene gene_name '
        '--sjdbOverhang 109 '
        '--outFileNamePrefix {params.outdir}/ '
        '&> {log}'

rule generate_gtf:
    input:
        os_genome = os_genome,
        os_gff = os_gff,
    output:
        os_gtf = os_gtf,
        cuffcomp_dir = temp('output/010_data/cuffcomp'),
        cuffcomp_gtf = temp('output/010_data/cuffcomp/'
                            'annot.cuffcomp.combined.gtf')
    threads:
        1
    log:
        'output/000_logs/010_prepare-data/generate_gtf.log'
    shell:
        'cuffcompare '
        '-s {input.os_genome} '
        '-CG '
        '-r {input.os_gff} '
        '-o {output.cuffcomp_dir}/annot.cuffcomp '
        '{input.os_gff} '
        '&> {log} ; '
        'sed \'/LOC_Os09g01000/d\' '
        '{output.cuffcomp_gtf} '
        ' | sed \'/LOC_Os09g00999/d\' '
        '> {output.os_gtf}'
