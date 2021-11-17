import pandas as pd

## functions ####
def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


## input data ###£
samplesfile = "samples.txt"
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)

rule all:
    input:
        "qc/multiqc_report.html",
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

rule align:
    input:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        gtf="genome/human.GRCh38.chr22.gtf",
        genome="genome/STARINDEX/Genome"
    output:
        "star/{sample}-{unit}.Aligned.out.sam",
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
            --quantMode GeneCounts \
            --outSAMtype SAM
        """

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
