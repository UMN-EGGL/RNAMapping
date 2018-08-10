
# export AWS_ACCESS_KEY=.....
# export AWS_SECRET_KEY=.....
# run the above lines in the shell before running snakemake 

# TO DO???
# -use config files? - could use https://github.com/UofABioinformaticsHub/RNAseq_snakemake
# -include log options for each rule?

# QUESTIONS
# -how many threads to use?
# -generate_index: GTF or GFF3? Overhang length for generating the splice
#   junction database?

import os
from snakemake.io import glob_wildcards 


from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider              
s3_key_id = os.environ.get('AWS_ACCESS_KEY')                                    
s3_access_key = os.environ.get('AWS_SECRET_KEY')                                
                                                                                
S3 = S3RemoteProvider(                                                          
    endpoint_url='https://s3.msi.umn.edu',                                      
    access_key_id=s3_key_id,                                                    
    secret_access_key=s3_access_key                                            
)



EQUCAB3_FASTA = S3.remote("FASTAs/EquCab2.nice.fna")
SAMPLES, = S3.glob_wildcards("fake_data/{sample}_R1_001.fastq")

#SAMPLES, *_ = GS.glob_wildcards(GS_PREFIX + '/griffithlab_brain_vs_uhr/HBR_UHR_ERCC_ds_10pc/{sample}.read1.fastq.gz')

rule all:   
    input:
        ...

rule trim_reads:
    input:
        fwd=expand("fake_data/{sample}_R1_001.fastq", sample=SAMPLES), 
        rev=expand("fake_data/{sample}_R2_002.fastq", sample=SAMPLES)
    output:
        fwd_trim="trimmed_data/{sample}_R1_001_trim.fastq", 
        rev_trim="trimmed_data/{sample}_R2_001_trim.fastq"
    message:
        "Executing AdapterRemoval on the following files {input}"
    shell:
        """
        AdapterRemoval \
#        --threads ??? \
        --file1 {input.fwd} \
        --file2 {input.rev} \
        --output1 {output.fwd_trim} \
        --output2 {output.rev_trim} \
        --trimns
        --trimqualities
        """
        
rule generate_index:
    input:
        genome = EQUCAB_FASTA,
#        gtf_file = ""
    output:
        "star_index"
    message:
        "Generating genome index for {input.genome} with {input.gtf_file}"
    shell:
        """
        STAR \
#        --runThreadN ??? \
        --runMode genomeGenerate \
        --genomeDir {output} \
        --genomeFastaFiles {input.genome} \
        --sjdbGTFfile ??? \
        --sjdbOverhang ??? \
        """
    
    
    
    
    
    

STAR  --runMode genomeGenerate --runThreadN <# cpus> --genomeDir <genome output directory> --genomeFastaFiles <input Genome FASTA file>

i.e. STAR  --runMode genomeGenerate --runThreadN 24 --genomeDir ./ --genomeFastaFiles mm9.fa
