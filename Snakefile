import pandas as pd

samplesfile = "samples.txt"



samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)

rule all:
    input:
        expand("trimmed/{samples.sample}-{samples.unit}.1.fastq",
            samples=samples.itertuples()),
        "genome/STARINDEX/Genome"



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

def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule cutadapt:
      input:
          get_fastq
      output:
          fastq1="trimmed/{sample}-{unit}.1.fastq",
          fastq2="trimmed/{sample}-{unit}.2.fastq",
          qc="trimmed/{sample}-{unit}.qc.txt"
      conda:
          "envs/trim.yaml"
      params:
          nodes=1,
          adapter="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
      log:
          "logs/cutadapt/{sample}-{unit}.log"
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
