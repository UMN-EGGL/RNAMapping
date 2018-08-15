import os


#EQUCAB3_FASTA = S3.remote("FASTAs/EquCab2.nice.fna")
#EQUCAB3_GTF = S3.remote("")
#FASTQ = S3.remote("RnaSeqData/Project_McCue_Project_022/*.fastq")
#SAMPLES, = S3.glob_wildcards(
#        "RnaSeqData/Project_McCue_Project_022/{id}_R1_001.fastq"
#)

EQUCAB3_FA = ""
EQUCAB3_GTF = ""
STAR_INDEX_DIR = ""
FASTQ_DIR = "real_data"
SAMPLES, = glob_wildcards(FASTQ_DIR + "/{id}_R1_001.fastq")
READS = ["R1", "R2"]
TRIMS = ["trim1", "trim2"]

# COMMENT: HAVE NOT MADE ANY CHANGES TO THREAD USAGE FOR ANY RULE

rule all:
    input:
        expand("qc/qc_raw/{id}_{read}_001_fastqc.html", id=SAMPLES, read=READS),
#        expand("trimmed_data/{id}_{trim}.fastq.gz", id=SAMPLES, trim=TRIMS),
        expand("qc/qc_trim/{id}_trim1_fastqc.html", id=SAMPLES),
	expand("qc/qc_trim/{id}_trim2_fastqc.html", id=SAMPLES)

rule trim_reads:
    input:
        R1 = FASTQ_DIR + "/{id}_R1_001.fastq",
        R2 = FASTQ_DIR + "/{id}_R2_001.fastq"
    output:
        R1 = "trimmed_data/{id}_trim1.fastq.gz",
        R2 = "trimmed_data/{id}_trim2.fastq.gz"
    message:
        "AdapterRemoval - removing adapters and low quality bases on {wildcards.id}"
    shell:
        """
        AdapterRemoval \
        --file1 {input.R1} \
        --file2 {input.R2} \
        --output1 {output.R1} \
        --output2 {output.R2} \
        --gzip \
        --trimns \
        --trimqualities \
        --minquality 10 \
        """


rule qc_trim:
    input:
        "trimmed_data/{id}.fastq.gz"
    params:
        out_dir = "qc/qc_trim"
    output:
        "qc/qc_trim/{id}_fastqc.html"
    message:
        "FastQC - performing quality control on trimmed {wildcards.id}"
    shell:
        """
        fastqc \
        -o {params.out_dir} \
        -f fastq \
        {input}
        """

rule qc_raw:
    input:
        FASTQ_DIR + "/{id}_{read}_001.fastq"
    params:
        out_dir = "qc/qc_raw"
    output:
        "qc/qc_raw/{id}_{read}_001_fastqc.html"
    message:
        "FastQC - performing quality control on {wildcards.id}_{wildcards.read}"
    shell:
        """
        fastqc \
        -o {params.out_dir} \
        -f fastq \
        {input}
        """

# QUESTIONS: DO WE NEED TO INCLUDE/DEAL WITH THE OVERHANG PARAMETER (--sjdbOverhang)?

rule STAR_index:
    input:
        genome_fa = EQUCAB3_FA,
        gtf = EQUCAB3_GTF
    output:
        STAR_INDEX_DIR
    message:
        "STAR - generating genome index for {input.genome_fa} with {input.gtf}"
    shell:
        """
        
        ""
