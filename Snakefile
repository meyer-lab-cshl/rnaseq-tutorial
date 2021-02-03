
import pandas as pd
samplesfile = "samples.txt"

samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)
print(samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna())

rule genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    threads: 1
    conda:
        "envs/align.yaml"
    params:
        saindex=11,
        overhang=75
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases {params.saindex} \
            --sjdbOverhang {params.overhang}
        """

def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


rule cutadapt:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        qc="trimmed/{sample}-{unit}.qc.txt"
    threads: 1
    conda:
        "envs/trim.yaml"
    params:
        adapter="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
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
