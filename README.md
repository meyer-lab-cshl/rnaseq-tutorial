# RNA-seq tutorial

We will use the workflow manager
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to develop
a RNAseq pipeline. One big part of developing successful pipelines is making
them reproducible and transferable, i.e. we should get the same results using the
pipeline on the same data, irrespective of which compute/system we analysed them
on. To achieve this, we will work with contained software environments, using
`mamba`/`conda`.

During the lectures, we will build our pipeline from scratch, starting with how
to write snakemake rules. In the end, we will have built an analysis pipeline for
RNAseq including alignment, quality control and differential expression analysis.
There is also a [homework](Homework.md) assignment that will ask you to extend
the pipeline, evaluate your results and learn how to generate an analysis report.

## Directory structure
The files you see in this directory are either part of the pipeline
([Snakefile](Snakefile),[scripts](scripts), [envs](envs)), or contain the example
data to run the analysis ([genome](genome), [reads](reads),
[samples.txt](samples.txt)).

## Set-up
To get started, please install the required software packages before the lectures.
The set-up is described below. If you run into issues with the installation,
please let me know BEFORE the first session, so we can resolve any problems
before the lectures.

(Don't worry if the setup does not make sense to you, or if you are not familiar
with the commands yet; we will cover it all in the lectures)

**1. Install a text editor**
Please make sure you have a text editor installed on your computer; if you do
not have one, give [VS Code](https://code.visualstudio.com/) a try!

**2. Install Miniforge**
If you do not have it installed already, we first need to install `miniforge` (Manufacturer's instruction here:
[miniforge](https://github.com/conda-forge/miniforge#mambaforge)).
For unix-like platforms (Mac and linux): open your terminal app (aka commad line, shell) and type:
```bash
curl -OL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
This will download and install miniforge. If all goes well you will see this in the end:

```bash
==> For changes to take effect, close and re-open your current shell. <==
```
Make sure you follow this advice and close and re-open your terminal. When you've done that, simply type `mamba` in your terminal
and you should see the `mamba` help manual, starting with

```bash
mamba                                                                                                                         [15:02:17]
usage: mamba [-h] [-V] command ...

conda is a tool for managing and deploying applications, environments and packages.
```

**2. Install Snakemake**
Following the recommodations by the
[snakemake developers](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html),
we will then set-up our snakemake environment:

```bash
mamba create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
```
mamba will check all packages and dependencies required to install snakemake and will
print these to the command line. 

At some point it will ask you `Confirm changes: [Y/n]`. Type `Y` and hit enter.

If it successfully finishes, you will see this on your screen:
```bash
Preparing transaction: done
Verifying transaction: done
Executing transaction: done

# To activate this environment, use
#
#     $ mamba activate snakemake-rnaseq
#
# To deactivate an active environment, use
#
#     $ mamba deactivate
```

**3. Activate the snakemake environment**
- type the following to activate the environment
```bash
mamba activate snakemake-rnaseq
```
- test the environment by typing
```bash
snakemake --version
```
You should see: `9.11.6`

- for those of you on a new Mac (using the Mx chips based on ARM64 architecture),
we also need to include this in our set-up to ensure all packages can be built
properly. You will only have to do this once:

```bash
conda config --env --set subdir osx-64
```

**4. Move into the rnaseq-tutorial directory**
- create a new folder for this tutorial. Do not include spaces in the folder
name/path,  e.g. do not call the folder ‘RNAseq tutorial’ use ‘rnaseq_tutorial’
instead.
- download and unzip the data and scripts from this repository by clicking on
the green ‘code’ button and selection ‘Download zip’. Alternatively, if you are
familiar with it, you can clone the repository.
- enter the directory using the command line `cd /path/2/rnaseq-tutorial`
  (where `path/2/` is the path to the directory where you saved the folder)


**5. Test the setup**
We currently are in a software environment that contains snakemake. We can
use this environment for any analysis that makes use of snakemake as a workflow
manager. That also means, we might use it for different pipelines that require
different software. To make everything transparent and reproducible, we
can tell snakemake itself that for each analysis type it should automatically
create a conda environment and handle the activation/deactivation, so we
do not have to worry about this. Here, we are testing this setup and the way
I have written it, all that snakemake will do here is create these environments
and use them to write some dummy text. If that works, we know we will not have
software related issues during the lecture and can instead focus on the
analysis, so let's do that!

- type: `snakemake --use-conda --cores 1 -s setup.smk`
- hit enter
- you should see something like this on your screen:
```
Building DAG of jobs...
Creating conda environment envs/pandas.yaml...
Downloading and installing remote packages.

...
```
(7 times for the different conda environments that snakemake will create)
followed by
```
Using shell: /bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Job stats:
job                count    min threads    max threads
---------------  -------  -------------  -------------
all                    1              1              1
count_matrix           1              1              1
cutadapt               1              1              1
deseq2                 1              1              1
generate_genome        1              1              1
index                  1              1              1
multiqc                1              1              1
rseqc_coverage         1              1              1
total                  8              1              1

```
- if all went well you should see this at the end:
```bash
Finished job 0.
8 of 8 steps (100%) done
```

If you do not see this, or get stuck anywhere before that, let me know and
we can have a look at it together - again, BEFORE class, so we can spend the
lecture time most efficiently.
