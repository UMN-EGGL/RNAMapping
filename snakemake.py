import os

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu', 
    access_key_id=s3_key_id, 
    secret_access_key=s3_access_key
)


SAMPLES, = S3.glob_wildcards('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq.gz')

# HAVE NOT MADE ANY CHANGES TO THREAD USAGE FOR ANY RULE

rule all:
    input:
        S3.remote(expand('qc/qc_raw/{sample}_fastqc.html', sample=SAMPLES))
        S3.remote(expand('qc/qc_trim/{sample}_fastqc.html', sample=SAMPLES))
        S3.remote(expand('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}.bam', sample=SAMPLES))

# DOES NOT DEAL WITH .discarded.gz, .settings, .signleton.truncated.gz

rule trim_reads:
    input:
        R1 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R1_001.fastq',
        R2 = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/fastq/{sample}_R2_001.fastq'
    output:
        R1 = S3.remote(temp('trimmed_data/{sample}_trim1.fastq.gz')),
        R2 = S3.remote(temp('trimmed_data/{sample}_trim2.fastq.gz'))
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

# DOES NOT STORE UNMAPPED READS (--outReadsUnmapped) nor log output

rule STAR_mapping:
    input:
        R1 = S3.remote('trimmed_data/{sample}_trim1.fastq.gz'),
        R2 = S3.remote('trimmed_data/{sample}_trim2.fastq.gz')
        star_index = S3.remote('HorseGeneAnnotation/public/refgen/GCF_002863925.1_EquCab3.0')
    params:
        out_prefix = S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}')
    output:
        S3.remote('HorseGeneAnnotation/private/sequence/RNASEQ/bam/{sample}.bam')
    message:
        'STAR - performing mapping on {wildcards.sample} with {input.star_index}'
    shell:
        '''
        STAR \
        --genomeDir {input.star_index} \
        --readFilesIn {input.R1} {input.R2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix {params.out_prefix} \
        --outSAMtype BAM Unsorted \
        '''
