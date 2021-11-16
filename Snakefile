rule all:
    input:
        "genome/STARINDEX/Genome",
        "trimmed/Id1_AA-rep1.1.fastq"



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
    shell:
          """
          STAR \
              --runMode genomeGenerate \
              --runThreadN 1 \
              --genomeFastaFiles {input.genome} \
              --sjdbGTFfile {input.gtf} \
              --genomeDir genome/STARINDEX \
              --genomeSAindexNbases {params.Nbases} \
              --sjdbOverhang {params.length}
          """

rule cutadapt:
      input:
          fastq1="reads/S01_S1_R1_001.fastq",
          fastq2="reads/S01_S1_R2_001.fastq",
      output:
          fastq1="trimmed/Id1_AA-rep1.1.fastq",
          fastq2="trimmed/Id1_AA-rep1.2.fastq",
          qc="trimmed/Id1_AA-rep1.qc.txt"
      conda:
          "envs/trim.yaml"
      params:
          nodes=1,
          adapter="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
      log:
          "logs/cutadapt/Id1_AA-rep1.log"
      shell:
          """
          cutadapt \
              -a {params.adapter} \
              -o {output.fastq1} \
              -p {output.fastq2} \
              -j {params.nodes} \
              {input} \
          > {output.qc}
          """
