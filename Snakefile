##### rules #####
rule generate_genome:
    input:
        fa="genome/human.GRCh38.chr22.fasta",
        gtf="genome/human.GRCh38.chr22.gtf"
    output:
        genome="genome/STARINDEX/Genome"
    conda:
        "envs/align.yaml"
    shell:
        """
        # Using STAR for genome generation
        STAR \
            --runMode genomeGenerate \
            --runThreadN 1 \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} \
            --genomeDir genome/STARINDEX \
            --genomeSAindexNbases 11 \
            --sjdbOverhang 75
        """
