import pandas as pd

## set up data ####
samplesfile = "samples.txt"
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)
#print(samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna())

## functions ####
def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

## rules ####

rule all:
    input:
        expand("star/{samples.sample}-{samples.unit}.Aligned.sortedByCoord.out.bam",
            samples=samples.itertuples()),
        "qc/multiqc_report.html"

rule generate_genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    conda:
        "envs/align.yaml"
    params:
        length=75,
        Nbases=11
    threads: 1
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases {params.Nbases} \
            --sjdbOverhang {params.length}
        """

rule cutadapt:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.R1.fastq",
        fastq2="trimmed/{sample}-{unit}.R2.fastq",
        qc="trimmed/{sample}-{unit}.qc.txt"
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    params:
        adapters="CTGACCTCAAGTCTGCACACGAGAAGGCTAG",
    threads: 1
    shell:
        """
        cutadapt \
            -a {params.adapters} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            -j {threads} \
            {input} \
        > {output.qc}
        """

rule multiqc:
    input:
        expand("trimmed/{samples.sample}-{samples.unit}.qc.txt",
            samples=samples.itertuples())
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    conda:
        "envs/multiqc.yaml"
    shell:
        """
        multiqc \
            --force \
            --export \
            --outdir qc \
            --filename multiqc_report.html \
            trimmed > {log}
        """

rule align:
    input:
        fastq1="trimmed/{sample}-{unit}.R1.fastq",
        fastq2="trimmed/{sample}-{unit}.R2.fastq",
        gtf="genome/human.GRCh38.chr22.gtf",
        genome="genome/STARINDEX/Genome"
    output:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
        "star/{sample}-{unit}.ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        indexdir="genome/STARINDEX"
    threads: 4
    conda:
        "envs/align.yaml"
    shell:
        """
        STAR \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.indexdir} \
            --readFilesIn {input.fastq1} {input.fastq2} \
            --outFileNamePrefix star/{wildcards.sample}-{wildcards.unit}. \
            --quantMode GeneCounts \
            --outSAMtype BAM SortedByCoordinate
        """
