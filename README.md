# RNA-seq tutorial

We will use the workflow manager [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to develop
a RNAseq pipeline. One big part of developing successful pipelines is making them reproducible and transferable,
i.e. we should get the same results using the pipeline on the same data, irrespective of which compute/system we
analysed them on. To achieve this, we will work with contained software environments, using `conda` (which you
should all have installed already when you installed [Anaconda](https://www.anaconda.com/products/individual).

To get started, please install the required software packages before the lectures.
The set-up is described below. If you run into issues with the installation, please let me know BEFORE the
first session, so we can resolve any problems before the lectures.

(Don't worry if the setup does not make sense to you, or if you are not familiar with the commands yet; we will
cover it all in the lectures)

## Set-up
**1. Install a text editor**
Please make sure you have a text editor installed on your computer; if you do
not have one, give [Atom}(https://atom.io/) a try!

**2. Create a conda environment for the tutorial.**
- Open your terminal app (aka commad line, shell)
- Type the following command:

```bash
conda create -c conda-forge -c bioconda -n snakemake-rnaseq snakemake
```
NB: this might take a while. Initially you should see:
```bash
Collecting package metadata (current_repodata.json): done
Solving environment: done
```
At some point it will ask you `Proceed ([y]/n)?`. Type `y` and hit enter.

If it successfully finishes, you will see this on your screen:
```bash
# To activate this environment, use
#
#     $ conda activate snakemake-rnaseq
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

**3. activate the snakemake environment** 
- type the following to activate the environment
```bash
conda activate snakemake-rnaseq
```
- test the environment by typing
```bash
snakemake --version
```
You should see: `5.28.0`

**4. Move into the rnaseq-tutorial directory**
- enter the directory using the command line `cd /path/2/rnaseq-tutorial`
  (where `path/2/` is the path to the directory where you saved the folder)

**5. Test the setup**
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


