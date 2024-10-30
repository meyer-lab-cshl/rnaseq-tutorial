import pandas as pd

samplesfile = "samples.txt"
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)
print(samples)

rule all:
    input:
        expand("trimmed/{samples.sample}-{samples.unit}.1.fastq",
            samples=samples.itertuples()),
        "qc/multiqc_report.html"  

##### rules #####
rule generate_genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    threads: 2
    conda:
        "envs/align.yaml"
    params:
        length=75,
        Nbases=11,
        genomedir="genome/STARINDEX"
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir {params.genomedir} \
            --genomeSAindexNbases {params.Nbases} \
            --sjdbOverhang {params.length}
        """
rule cutadapt:
    input:
        fastq1="reads/{sample}-{unit}.R1.fastq",
        fastq2="reads/{sample}-{unit}.R2.fastq",
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        qc="trimmed/{sample}-{unit}.qc.txt"
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    shell:
        """
        cutadapt \
            -a CTGACCTCAAGTCTGCACACGAGAAGGCTAG \
            -o {output.fastq1} \
            -p {output.fastq2} \
            -j 1 \
            {input} \
        > {output.qc}
        """

#################################################
# multi qc to visualise results of trim
#################################################
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