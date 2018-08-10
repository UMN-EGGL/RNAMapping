import os


#EQUCAB3_FASTA = S3.remote("FASTAs/EquCab2.nice.fna")
#EQUCAB3_GTF = S3.remote("")
#FASTQ = S3.remote("RnaSeqData/Project_McCue_Project_022/*.fastq")
#SAMPLES, = S3.glob_wildcards(
#        "RnaSeqData/Project_McCue_Project_022/{id}_R1_001.fastq"
#)

FASTQ_DIR = "real_data"
SAMPLES, = glob_wildcards(FASTQ_DIR + "/{id}_R1_001.fastq")
READS = ["R1", "R2"]
TRIMS = ["trim1", "trim2"]

rule all:
    input:
        expand("qc/qc_raw/{id}_{read}_001_fastqc.html", id=SAMPLES, read=READS),
        expand("trimmed_data/{id}/{id}_R1_001_trim1.fastq.gz", id=SAMPLES),
        expand("qc/qc_trim/{id}_{read}_001_{trim}_fasqc.html", id=SAMPLES, read=READS, trim=TRIMS)     

rule trim_reads:
    input:
        R1 = FASTQ_DIR + "/{id}_R1_001.fastq",
        R2 = FASTQ_DIR + "/{id}_R2_001.fastq"
    output:
        R1 = "trimmed_data/{id}/{id}_R1_001_trim1.fastq.gz",
        R2 = "trimmed_data/{id}/{id}_R2_001_trim2.fastq.gz"
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
        "trimmed_data/{id}/{id}_{read}_001_{trim}.fastq.gz"
    params:
        out_dir = "qc/qc_trim"
    output:
        "qc/qc_trim/{id}_{read}_001_{trim}_fastqc.html"
    shell:
        """
        fastqc \
        -o {params.out_dir} \
        -f fastq \
        {input}
        """

rule qc_raw:
    input:
        FASTQ_DIR + "/{id}_{read}_001.fastq",
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


