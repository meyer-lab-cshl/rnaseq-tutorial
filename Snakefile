import pandas as pd

def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


# parameters
samplesfile = "samples.txt"

# read samples to analyse
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)
print(samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna())


##### target rule
rule all:
    input:
        expand("trimmed/{samples.sample}-{samples.unit}.1.fastq",
            samples=samples.itertuples())

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
