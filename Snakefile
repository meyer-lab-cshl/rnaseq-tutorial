
##### rules #####
rule generate_genome:
    input:
        genome="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        "genome/STARINDEX/Genome"
    conda:
        "envs/align.yaml"
    shell:
        """
          STAR \
              --runMode genomeGenerate \
              --runThreadN 1 \
              --genomeFastaFiles {input.genome} \
              --sjdbGTFfile {input.gtf} \
              --genomeDir genome/STARINDEX \
              --genomeSAindexNbases 11 \
              --sjdbOverhang 75
        """


rule cutadapt:
    input:
        fastq1="reads/Id1_AA-rep1.R1.fastq",
        fastq2="reads/Id1_AA-rep1.R2.fastq",
    output:
        fastq1="trimmed/Id1_AA-rep1.R1.fastq",
        fastq2="trimmed/Id1_AA-rep1.R2.fastq",
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
