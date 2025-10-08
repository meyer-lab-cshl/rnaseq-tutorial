import pandas as pd

##### data ####
DESIGN="~ condition"
SPECIES="human"
SAMPLESFILE="samples.txt"

samples = pd.read_table(SAMPLESFILE).set_index(["sample", "unit"], drop=False)
print(samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna())


##### functions ####
def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


##### rules #####
rule all:
    input:
        expand("star/{samples.sample}-{samples.unit}.Aligned.sortedByCoord.out.bam",
            samples=samples.itertuples()),
        expand(["results/diffexp/{contrast}.diffexp.txt",
                "results/diffexp/{contrast}.ma-plot.pdf"],
               contrast=['AA', 'control']),
        "qc/multiqc_report.html"



rule generate_genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    threads: 1
    conda:
        "envs/align.yaml"
    params:
        length=75,
        Nbases=11
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
    params:
        adapters="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
    threads: 1
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
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
            --sjdbGTFfile {input.gtf} \
            --quantMode GeneCounts \
            --outSAMtype BAM SortedByCoordinate
        """

rule count_matrix:
    input:
        expand("star/{samples.sample}-{samples.unit}.ReadsPerGene.out.tab",
            samples=samples.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=samples['sample'].tolist(),
        strand="reverse"
    log:
        "logs/counts/count_matrix.log"
    conda:
       "envs/pandas.yaml"
    script:
        "scripts/count-matrix.py"


rule setup_de:
    input:
        counts="counts/all.tsv",
        annotation="genome/ENSEMBL_GRCh38p13.txt",
        samples=SAMPLESFILE
    output:
        dds="deseq2/all.rds"
    params:
        species=SPECIES,
        design=DESIGN,
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/setup.log"
    script:
        "scripts/setup_deseq2.R"


rule deseq2:
    input:
        dds="deseq2/all.rds",
    output:
        table="results/diffexp/{contrast}.diffexp.txt",
        ma_plot="results/diffexp/{contrast}.ma-plot.pdf",
        up="results/diffexp/deg-sig-up_{contrast}.csv",
        down="results/diffexp/deg-sig-down_{contrast}.csv"
    params:
        contrast=['AA', 'control'],
        design=DESIGN,
        samples=SAMPLESFILE
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/deseq2.R"

rule multiqc:
    input:
        expand("star/{samples.sample}-{samples.unit}.Aligned.sortedByCoord.out.bam",
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
            trimmed star > {log}
        """
