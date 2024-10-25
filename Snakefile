#################################################
# Skeleton of a rule
#################################################
rule all:
    input:
        "trimmed/Id1_AA-rep1.1.fastq",

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
        fastq1="reads/Id1_AA-rep1.R1.fastq",
        fastq2="reads/Id1_AA-rep1.R2.fastq",
    output:
        fastq1="trimmed/Id1_AA-rep1.1.fastq",
        fastq2="trimmed/Id1_AA-rep1.2.fastq",
        qc="trimmed/Id1_AA-rep1.qc.txt"
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/Id1_AA-rep1.log"
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