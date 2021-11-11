##### rules #####
rule all:
    input:
        "results/align.txt",
        "results/trim.txt",
        "results/index.txt",
        "results/rseqc.txt",
        "results/pandas.txt",
        "results/deseq2.txt",
        "results/multiqc.txt",


rule generate_genome:
    conda:
        "envs/align.yaml"
    output:
        "results/align.txt",
    shell:
        """
        star --version > results/align.txt
        echo 'align environment works' >> results/align.txt
        """

rule cutadapt:
    conda:
        "envs/trim.yaml"
    output:
        "results/trim.txt",
    shell:
        """
        cutadapt --version > results/trim.txt
        echo 'trim environment works' >> results/trim.txt
        """

rule index:
    conda:
        "envs/index.yaml"
    output:
        "results/index.txt",
    shell:
        """
        samtools --version > results/index.txt
        echo 'index environment works' >> results/index.txt
        """

rule rseqc_coverage:
    conda:
        "envs/rseqc.yaml"
    output:
        "results/rseqc.txt",
    shell:
        """
        geneBody_coverage.py --version > results/rseqc.txt
        infer_experiment.py --version >> results/rseqc.txt
        echo 'rseqc environment works' >> results/rseqc.txt
        """

rule count_matrix:
    conda:
       "envs/pandas.yaml"
    output:
        "results/pandas.txt",
    shell:
        """
        python -c 'import pandas as pd' > results/pandas.txt
        echo 'pandas environment works' >> results/pandas.txt
        """

rule deseq2:
    conda:
        "envs/deseq2.yaml"
    output:
        "results/deseq2.txt",
    shell:
        """
        Rscript -e \
        'library(tidyverse); library(DESeq2); library(biomaRt); library(cowplot)' \
        > results/deseq2.txt
        echo 'deseq2 environment works' >> results/deseq2.txt
        """

rule multiqc:
    conda:
        "envs/multiqc.yaml"
    output:
        "results/multiqc.txt",
    shell:
        """
        multiqc --version > results/multiqc.txt
        echo 'multiqc environment works' >> results/multiqc.txt
        """



