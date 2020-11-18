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
        fastq1="reads/S01_S1_R1_001.fastq",
        fastq2="reads/S01_S1_R2_001.fastq",
    output:
        fastq1="trimmed/Id1_AA-rep1.1.fastq",
        fastq2="trimmed/Id1_AA-rep1.2.fastq",
        qc="trimmed/Id1_AA-rep1.qc.txt"
    params:
        adapter="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
    threads: 1
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/Id1_AA-rep1.log"
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
