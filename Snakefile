import pandas as pd

def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


# parameters
samplesfile = "samples.txt"

# read samples to analyse
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)

##### target rule
rule all:
    input:
        expand("star/{samples.sample}-{samples.unit}.Aligned.out.sam",
            samples=samples.itertuples()),
        "qc/multiqc_report.html"

##### rules #####
rule generate_genome:
    input:
        fa="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        genome="genome/STARINDEX/Genome"
    threads: 1
    params:
        saindex = 11,
        overhang = 75
    conda:
        "envs/align.yaml"
    shell:
        """
        # Using STAR for genome generation
        STAR \
            --runMode genomeGenerate \
            --runThreadN  {threads} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases {params.saindex} \
            --sjdbOverhang {params.overhang}
        """

rule cutadapt:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        adapter="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
    threads: 1
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt \
            -a {params.adapter} \
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
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        genome="genome/STARINDEX/Genome"
    output:
        "star/{sample}-{unit}.Aligned.out.sam",
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
            --outSAMtype SAM
        """
