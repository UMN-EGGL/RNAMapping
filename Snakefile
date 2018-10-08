import os
import glob
import re
from collections import Counter


from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu', 
    access_key_id=s3_key_id, 
    secret_access_key=s3_access_key
)

se_samples_tmp = []
for path in [x for x in S3._s3c.list_keys('HorseGeneAnnotation') if 'fastq' in x]: 
    dir_path, se_file = os.path.split(path)
    se_file = re.search(r'^(.*?)(_R[1-2]_)',se_file).group(1)
    se_samples_tmp.append(se_file)
    
SE_SAMPLES = [k for k, v in Counter(se_samples_tmp).items() if v == 1]

SAMPLES, = S3.glob_wildcards('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R2_001.fastq.gz')
configfile: "config.yaml"

rule all:
    input:
        #S3.remote(expand('qc/qc_raw/{sample}_fastqc.html', sample=SAMPLES)),
        #S3.remote(expand('qc/qc_trim/{sample}_fastqc.html', sample=SAMPLES)),
        #S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_Aligned.out.bam', sample=SAMPLES)),
        S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_se_Aligned.out.bam', sample=SE_SAMPLES)),
        #gff = S3.remote( expand("HorseGeneAnnotation/public/refgen/{GCF}/GFF/{sample}.gff" ,sample=SAMPLES,GCF=config['GCF']))

# ----------------------------------------------------------
#       Trimming
# ----------------------------------------------------------

rule trim_reads:
    input:
        R1 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz'),
        R2 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R2_001.fastq.gz')
    output:
        R1 = temp('trimmed_data/{sample}_trim1.fastq.gz'),
        R2 = temp('trimmed_data/{sample}_trim2.fastq.gz')
    message:
        'AdapterRemoval - removing adapters and low quality bases on {wildcards.sample}'
    shell:
        '''
        AdapterRemoval \
        --file1 {input.R1} \
        --file2 {input.R2} \
        --output1 {output.R1} \
        --output2 {output.R2} \
        --gzip \
        --trimns \
        --trimqualities \
        --minquality 10 \
        '''

rule get_se_fastqs:
    input:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz')
    output:
        'local/HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz'
    run:
        shell('cp {input[0]} {output[0]}')

rule trim_se_read:
    input:
        R1 = 'local/HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz' 
        #R1 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz')
    output:
        R1 = temp('trimmed_data/{sample}_se_trim.fastq.gz')
    message:
        'AdapterRemoval - removing adapters and low quality bases on SE reads {wildcards.sample}'
    shell:
        '''
        AdapterRemoval \
        --file1 {input.R1} \
        --output1 {output.R1} \
        --gzip \
        --trimns \
        --trimqualities \
        --minquality 10 \
        '''

# ----------------------------------------------------------
#       QC
# ----------------------------------------------------------

rule qc_trim:
    input:
        S3.remote('trimmed_data/{sample}.fastq.gz')
    params:
        out_dir = S3.remote('qc/qc_trim')
    output:
        S3.remote('qc/qc_trim/{sample}_fastqc.html')
    message:
        'FastQC - performing quality control on trimmed {wildcards.sample}'
    shell:
        '''
        fastqc \
        -o {params.out_dir} \
        -f fastq \
        {input}
        '''

rule qc_raw:
    input:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}.fastq.gz')
    params:
        out_dir = S3.remote('qc/qc_raw')
    output:
        S3.remote('qc/qc_raw/{sample}_fastqc.html')
    message:
        'FastQC - performing quality control on {wildcards.sample}'
    shell:
        '''
        fastqc \
        -o {params.out_dir} \
        -f fastq \
        {input}
        '''

# ----------------------------------------------------------
#       STAR Mapping
# ----------------------------------------------------------

rule download_STAR:
    input:
        expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done',GCF=config['GCF'])

#DOES NOT STORE UNMAPPED READS (--outReadsUnmapped) nor log output
rule STAR_index:
    input:
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/Genome',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/SA',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/SAindex',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrLength.txt',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrName.txt',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrNameLength.txt',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/chrStart.txt',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/exonGeTrInfo.tab',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/exonInfo.tab',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/geneInfo.tab',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/genomeParameters.txt',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbInfo.txt',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbList.fromGTF.out.tab',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/sjdbList.out.tab',GCF=config['GCF']),keep_local=True),
        S3.remote(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/transcriptInfo.tab',GCF=config['GCF']),keep_local=True) 
    output:
        touch(expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done',GCF=config['GCF']))

rule STAR_mapping:
    input:
        R1 = 'trimmed_data/{sample}_trim1.fastq.gz',
        R2 = 'trimmed_data/{sample}_trim2.fastq.gz',
        star_index = expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES/download.done',GCF=config['GCF'])
    output:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_Aligned.out.bam')
    params:
        out_prefix = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_'),
        star_index = expand('HorseGeneAnnotation/public/refgen/{GCF}/STAR_INDICES',GCF=config['GCF'])
    message:
        'STAR - Creating: {output} '
    run:
        assert os.path.exists(input.star_index)
        shell('''
        STAR \
        --genomeDir {params.star_index} \
        --genomeLoad LoadAndKeep \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.out_prefix} \
        --outSAMtype BAM Unsorted \
        ''')

rule STAR_mapping_se:
    input:
        R1 = 'trimmed_data/{sample}_se_trim.fastq.gz',
        star_index = 'HorseGeneAnnotation/public/refgen/GCF_002863925.1_EquCab3.0/STAR_INDICES/download.done'
    params:
        out_prefix = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_se_'),
        star_index = 'HorseGeneAnnotation/public/refgen/GCF_002863925.1_EquCab3.0/STAR_INDICES'
    output:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_se_Aligned.out.bam')
    message:
        'STAR - Creating: {output} '
    run:
        assert os.path.exists(input.star_index)
        shell('''
        STAR \
        --genomeDir {params.star_index} \
        --genomeLoad LoadAndKeep \
        --readFilesIn {input.R1} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.out_prefix} \
        --outSAMtype BAM Unsorted \
        ''')

# ----------------------------------------------------------
#       StringTie
# ----------------------------------------------------------

# 90M_ATCACG_L004

rule run_stringtie:
    input:
        bam = S3.remote(ancient('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}_Aligned.out.bam')),
        gff = S3.remote(expand(ancient("HorseGeneAnnotation/public/refgen/{GCF}/{GCF}_genomic.nice.gff.gz"), GCF=config['GCF']))
    output:
        gff = S3.remote( f"HorseGeneAnnotation/public/refgen/{config['GCF']}/GFF/{{sample}}.gff" )
    run:
        shell('''
        stringtie \
        {input.bam} \
        -G {input.gff} \
        -o {output.gff}
        ''')

