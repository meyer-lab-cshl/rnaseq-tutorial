import pandas as pd

## functions ####
def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


## input data ###
samplesfile = "samples.txt"
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)

rule all:
    input:
        "qc/multiqc_report.html",
        "results/diffexp/AA_vs_control.diffexp.txt"



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
        samples="samples.txt",
        annotation="genome/ENSEMBL_GRCh38p13.txt"
    output:
        dds="deseq2/all.rds"
    params:
        species="human",
        design="~ condition",
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
        design="~ condition",
        samples="samples.txt"
    conda:
        "envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    script:
        "scripts/deseq2.R"

rule index:
    input:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
    conda:
        "envs/index.yaml"
    shell:
        "samtools index {input}"


rule rseqc_coverage:
    input:
        bed="genome/human.GRCh38.chr22.bed",
        bam="star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
        bai="star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "qc/rseqc/{sample}-{unit}.geneBodyCoverage.txt"
    log:
        "logs/rseqc/rseqc_coverage/{sample}-{unit}.log"
    conda:
        "envs/rseqc.yaml"
    shell:
        """
        geneBody_coverage.py \
            -r {input.bed} \
            -i {input.bam} \
            -o qc/rseqc/{wildcards.sample}-{wildcards.unit} 2> {log}
        """

rule multiqc:
    input:
        expand("qc/rseqc/{samples.sample}-{samples.unit}.geneBodyCoverage.txt",
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
            trimmed star qc/rseqc > {log}
        """
