# RNA-seq tutorial

We will use the workflow manager [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to develop
a RNAseq pipeline. To get started, please install the required software packages before the lectures.
The set-up is described below. If you run into issues with the installation, please let me know BEFORE the
first session, so we can resolve any problems before the lectures. 

## Set-up
**1. Create a conda environment for the tutorial.**
- Open your terminal app (aka commad line, shell)
- Type the following command:

```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake-rnaseq
```
- activate the snakemake environment with 
```bash
conda activate snakemake-rnaseq
```

**2. Save rnaseq-tutorial.zip**
- unzip
- move into the directory using the command line `cd /path/2/rnaseq-tutorial`
  (where `path/2/` is the path to the directory where you saved the folder)

**3. Test the setup**
- type: `snakemake --use-conda --cores 1`
- hit enter
- you should see something like this on your screen:
```
Building DAG of jobs...
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job counts:
	count	jobs
	1	all
	1	count_matrix
	1	cutadapt
	1	deseq2
	1	generate_genome
	1	index
	1	multiqc
	1	rseqc_coverage
	8
```
- if all went well you should see this at the end:
```bash
Finished job 0.
8 of 8 steps (100%) done
```

NB: You can ignore ERRORs like this:
```
ERROR: This cross-compiler package contains no program /Users/hannah/teaching/rnaseq/.snakemake/conda/381cf40c/bin/x86_64-apple-darwin13.4.0-ar
ERROR: activate_clang_osx-64.sh failed, see above for details
ERROR: This cross-compiler package contains no program /Users/hannah/teaching/rnaseq/.snakemake/conda/381cf40c/bin/x86_64-apple-darwin13.4.0-clang++
ERROR: activate_clangxx_osx-64.sh failed, see above for details
```

## Analysis

## Results
 
