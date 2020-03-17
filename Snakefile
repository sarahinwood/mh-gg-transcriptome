#!/usr/bin/env python3

import pathlib2
import pandas
import os
import snap

######trinity simg download - add to script and try run

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def find_read_files(read_dir):
#Make list of files
    path_generator = os.walk(read_dir, followlinks = True)
    my_files = list((dirpath, filenames)
        for (dirpath, dirname, filenames)
        in path_generator)
#Make new dictionary & populate with files (flowcell = key)
    my_fastq_files = {}
    for dirpath, filenames in my_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz'):
                my_flowcell = pathlib2.Path(dirpath).name
                my_fastq = str(pathlib2.Path(dirpath,filename))
                if my_flowcell in my_fastq_files:
                    my_fastq_files[my_flowcell].append(my_fastq)
                else:
                    my_fastq_files[my_flowcell]= []
                    my_fastq_files[my_flowcell].append(my_fastq)
    return(my_fastq_files)

def sample_name_to_fastq(wildcards):
    sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
    sample_id = sample_row.iloc[-1]['OGF_sample_ID']
    sample_flowcell = sample_row.iloc[-1]['Flow_cell']
    sample_all_fastq = [x for x in all_fastq[sample_flowcell]
                        if '-{}-'.format(sample_id) in x]
    sample_r1 = sorted(list(x for x in sample_all_fastq
                            if '_R1_' in os.path.basename(x)))
    sample_r2 = sorted(list(x for x in sample_all_fastq
                            if '_R2_' in os.path.basename(x)))
    return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

read_dir = 'data/reads'

sample_key_file = 'data/sample_key.csv'

star_reference_folder = 'output/star/star_reference'

#containers
busco_container = 'shub://TomHarrop/singularity-containers:busco_3.0.2'
tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
star_container = 'shub://TomHarrop/singularity-containers:star_2.7.0c'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'

#########
# SETUP #
#########

# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

sample_key = pandas.read_csv(sample_key_file)
all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
    input:
        'output/trinity/Trinity-GG.fasta',
        expand('output/busco/run_{filter}/full_table_{filter}.tsv',
               filter=['expression', 'length']),
        'output/trinity_stats/stats.txt',
        'output/trinity_stats/xn50.out.txt',
        'output/trinity_stats/bowtie2_alignment_stats.txt',
        'output/transrate/Trinity-GG/contigs.csv',
        'output/trinotate/trinotate/Trinotate.sqlite'

rule trinotate:
    input:
        fasta = 'output/trinity/Trinity-GG.fasta',
        blastdb = 'bin/trinotate/db/uniprot_sprot.pep',
        hmmerdb = 'bin/trinotate/db/Pfam-A.hmm',
        sqldb = 'bin/trinotate/db/Trinotate.sqlite'
    output:
        'output/trinotate/trinotate/trinotate_annotation_report.txt',
        'output/trinotate/trinotate/Trinotate.sqlite'
    params:
        wd = 'output/trinotate'
    threads:
        20
    log:
        'output/logs/trinotate.log'
    shell:
        'trinotate_pipeline '
        '--trinity_fasta {input.fasta} '
        '--blast_db {input.blastdb} '
        '--hmmer_db {input.hmmerdb} '
        '--sqlite_db {input.sqldb} '
        '--outdir {params.wd} '
        '--threads {threads} '
        '2> {log}'

rule busco:
    input:
        filtered_fasta = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.fasta',
        lineage = 'data/hymenoptera_odb9'
    output:
        'output/busco/run_{filter}/full_table_{filter}.tsv'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                            'busco_{filter}.log'))
    params:
        wd = 'output/busco',
        filtered_fasta = lambda wildcards, input: resolve_path(input.filtered_fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage)
    threads:
        20
    singularity:
        busco_container
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--in {params.filtered_fasta} '
        '--out {wildcards.filter} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--species nasonia '
        '--mode transcriptome '
        '-f '
        '2> {log} '

rule transrate:
    input:
        transcriptome = 'output/trinity/Trinity-GG.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples)
    output:
        'output/transrate/Trinity-GG/contigs.csv'
    log:
        'output/logs/transrate.log'
    params:
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right))),
        outdir = 'output/transrate/'
    threads:
        50
    shell:
        'bin/transrate/transrate '
        '--assembly {input.transcriptome} '
        '--left {params.left} '
        '--right {params.right} '
        '--output {params.outdir} '
        '--threads {threads} '
        '--loglevel error '
        '2> {log}'

rule bowtie2_alignment_stats:
    input:
        transcriptome = 'output/trinity/Trinity-GG.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples)
    output:
        alignment_stats = 'output/trinity_stats/bowtie2_alignment_stats.txt'
    params:
        index_basename = 'output/trinity_stats/Trinity.fasta.index',
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right)))
    threads:
        50
    shell:
        'bowtie2-build '
        '{input.transcriptome} '
        '{params.index_basename} || exit 1 ; '
        'bowtie2 '
        '-p 10 '
        '-q '
        '--threads {threads} '
        '-x {params.index_basename} '
        '-1 {params.left} '
        '-2 {params.right} '
        '1> /dev/null 2> {output.alignment_stats}'

rule filter_trinity_isoforms:
    input:
        transcriptome = 'output/trinity/Trinity-GG.fasta',
        isoforms = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.txt'
    output:
        sorted_fasta = 'output/trinity_filtered_isoforms/isoforms_by_{filter}.fasta'
    log:
        'output/logs/filter_isoforms_by_{filter}.log'
    singularity:
        bbduk_container
    shell:
        'filterbyname.sh '
        'in={input.transcriptome} '
        'include=t '
        'names={input.isoforms} '
        'out={output.sorted_fasta} ' 
        '2> {log}'

rule sort_isoforms_r:
    input:
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        expression = 'output/trinity_filtered_isoforms/isoforms_by_expression.txt',
        length = 'output/trinity_filtered_isoforms/isoforms_by_length.txt'
    singularity:
        tidyverse_container
    log:
        'output/logs/sort_isoforms_r.log'
    script:
        'scripts/sort_isoforms_mh.R'

rule ExN50_stats:
    input:
        abundance = 'output/trinity_abundance/RSEM.isoform.TPM.not_cross_norm',
        transcriptome = 'output/trinity/Trinity-GG.fasta'
    output:
        ExN50_stats = 'output/trinity_stats/xn50.out.txt'
    log:
        'output/logs/xn50.err.txt'
    shell:
        'singularity exec -e singularity/trinityrnaseq.v2.9.1.simg /usr/local/bin/trinityrnaseq/util/misc/contig_ExN50_statistic.pl '
        '{input.abundance} '
        '{input.transcriptome} '
        '>{output.ExN50_stats} '
        '2>{log}'

rule trinity_stats:
    input:
        transcriptome = 'output/trinity/Trinity-GG.fasta'
    output:
        stats = 'output/trinity_stats/stats.txt'
    log:
        'output/logs/trinity_stats.log'
    shell:
        'singularity exec -e singularity/trinityrnaseq.v2.9.1.simg /usr/local/bin/trinityrnaseq/util/TrinityStats.pl '
        '{input.transcriptome} '
        '>{output.stats} '
        '2>{log}'

rule trinity_abundance_to_matrix:
    input:
        gt_map = 'output/trinity/Trinity-GG.fasta.gene_trans_map',
        abundance = 'output/trinity_abundance/RSEM.isoforms.results'
    output:
        'output/trinity_abundance/RSEM.isoform.counts.matrix',
        'output/trinity_abundance/RSEM.isoform.TPM.not_cross_norm'
    params:
        prefix = 'output/trinity_abundance/RSEM'
    log:
        'output/logs/abundance_estimates_to_matrix.log'
    shell:
        'singularity exec -e singularity/trinityrnaseq.v2.9.1.simg /usr/local/bin/trinityrnaseq/util/abundance_estimates_to_matrix.pl '
        '--est_method RSEM '
        '--cross_sample_norm none '
        '--out_prefix {params.prefix} '
        '--gene_trans_map {input.gt_map} '
        '{input.abundance} '
        '2> {log}'

rule trinity_abundance:
    input:
        transcripts = 'output/trinity/Trinity-GG.fasta',
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples)
    output:
        'output/trinity_abundance/RSEM.isoforms.results'
    threads:
        20
    log:
        'output/logs/trinity_abundance.log'
    params:
        outdir = 'output/trinity_abundance',
        left = lambda wildcards, input: ','.join(sorted(set(input.left))),
        right = lambda wildcards, input: ','.join(sorted(set(input.right)))
    shell:
        'singularity exec -e singularity/trinityrnaseq.v2.9.1.simg /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl '
        '--transcripts {input.transcripts} '
        '--seqType fq '
        '--est_method RSEM '
        '--aln_method bowtie2 '
        '--output_dir {params.outdir} '
        '--prep_reference '
        '--SS_lib_type RF '
        '--thread_count {threads} '
        '--trinity_mode '
        '--left {params.left} '
        '--right {params.right} '
        '2> {log}'

rule trinity_genome_guided:
    input:
        bam = 'output/star/star_pass2/Aligned.sortedByCoord.out.bam'
    output:
        fasta = 'output/trinity/Trinity-GG.fasta',
        gene_to_transcript = 'output/trinity/Trinity-GG.fasta.gene_trans_map'
    params:
        outdir = 'output/trinity'
    threads:
        30
    log:
        'output/logs/trinity.log'
    shell:
        'singularity exec -e singularity/trinityrnaseq.v2.9.1.simg Trinity '
        '--genome_guided_bam {input.bam} '
        '--genome_guided_max_intron 10000 '
        '--max_memory 300G '
        '--CPU {threads} '
        '--SS_lib_type RF '
        '--output {params.outdir} '
        '--seqType fq '
        '2> {log}'

##map reads to mh hic genome
rule star_second_pass:
    input:
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples),
        junctions = 'output/star/star_pass1/SJ.out.tab'
    output:
        bam = 'output/star/star_pass2/Aligned.sortedByCoord.out.bam'
    threads:
        30
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star/star_pass2/'
    log:
        'output/logs/star/star_pass2.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outBAMcompression 10 '
        '--readFilesIn output/bbduk_trim/Mhyp_abdo_r1.fq.gz,output/bbduk_trim/Mhyp_head_r1.fq.gz,output/bbduk_trim/Mhyp_ovary1_r1.fq.gz,output/bbduk_trim/Mhyp_ovary2_r1.fq.gz,output/bbduk_trim/Mhyp_sting_r1.fq.gz output/bbduk_trim/Mhyp_abdo_r2.fq.gz,output/bbduk_trim/Mhyp_head_r2.fq.gz,output/bbduk_trim/Mhyp_ovary1_r2.fq.gz,output/bbduk_trim/Mhyp_ovary2_r2.fq.gz,output/bbduk_trim/Mhyp_sting_r2.fq.gz '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '2> {log}'

##input being expanded doesn't pass list as needed
rule star_first_pass:
    input:
        left = expand('output/bbduk_trim/{sample}_r1.fq.gz', sample=all_samples),
        right = expand('output/bbduk_trim/{sample}_r2.fq.gz', sample=all_samples),
        star_reference = 'output/star/star_reference/Genome'
    output:
        sjdb = 'output/star/star_pass1/SJ.out.tab'
    params:
        genome_dir = star_reference_folder,
        prefix = 'output/star/star_pass1/'
    threads:
        30
    log:
        'output/logs/star/star_pass1.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '
        '--readFilesIn output/bbduk_trim/Mhyp_abdo_r1.fq.gz,output/bbduk_trim/Mhyp_head_r1.fq.gz,output/bbduk_trim/Mhyp_ovary1_r1.fq.gz,output/bbduk_trim/Mhyp_ovary2_r1.fq.gz,output/bbduk_trim/Mhyp_sting_r1.fq.gz output/bbduk_trim/Mhyp_abdo_r2.fq.gz,output/bbduk_trim/Mhyp_head_r2.fq.gz,output/bbduk_trim/Mhyp_ovary1_r2.fq.gz,output/bbduk_trim/Mhyp_ovary2_r2.fq.gz,output/bbduk_trim/Mhyp_sting_r2.fq.gz '
        '--readFilesCommand zcat '
        '--outFileNamePrefix {params.prefix} '
        '2> {log}'

rule star_reference:
    input:
        mh_genome = 'data/Mh_Hi-C_PGA_assembly.fasta'
    output:
        'output/star/star_reference/Genome'
    params:
        genome_dir = star_reference_folder
    threads:
        20
    log:
        'output/logs/star/star_reference.log'
    singularity:
        star_container
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.genome_dir} '
        '--genomeFastaFiles {input.mh_genome} '
        '2> {log} '

