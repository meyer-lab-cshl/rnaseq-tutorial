import pandas as pd
samplesfile = "samples.txt"
samples = pd.read_table(samplesfile).set_index(["sample", "unit"], drop=False)
print(samples.loc[('Id1_AA', 'rep1'), ["fq1", "fq2"]].dropna())

def get_fastq(wildcards):
    return samples.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule all:
    input:
        "genome/STARINDEX/Genome",
        "qc/multiqc_report.html",
        #expand("star/{samples.sample}-{samples.unit}.Aligned.sortedByCoord.out.bam",
        #    samples=samples.itertuples())
        expand("qc/rseqc/{samples.sample}-{samples.unit}.geneBodyCoverage.txt",
        samples=samples.itertuples())



rule build_genome:
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
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq",
        fastq2="trimmed/{sample}-{unit}.2.fastq",
        qc="trimmed/{sample}-{unit}.qc.txt"
    threads: 1
    params:
        adapter="CTGACCTCAAGTCTGCACACGAGAAGGCTAG"
    conda:
        "envs/trim.yaml"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
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

#################################################
# multi qc to visualise results of trim
#################################################
rule multiqc:
    input:
        expand("trimmed/{samples.sample}-{samples.unit}.qc.txt",
            samples=samples.itertuples()),
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
            trimmed star> {log}
        """

#################################################
# Align; BAM as output; one file only
#################################################

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


#################################################
# Helper rule for downstream commands: index bam
#################################################
rule index:
    input:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam",
    output:
        "star/{sample}-{unit}.Aligned.sortedByCoord.out.bam.bai",
    conda:
        "envs/index.yaml"
    shell:
        "samtools index {input}"

#################################################
# More QC: Gene body coverage
#################################################
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
