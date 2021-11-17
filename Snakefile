import pandas as pd

samplesfile = "samples.txt"



samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)
print(samples)
print(samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna())



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
          samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna()
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
