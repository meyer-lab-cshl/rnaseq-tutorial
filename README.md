# RNA-seq tutorial

We will use the workflow manager
[snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to develop
a RNAseq pipeline. One big part of developing successful pipelines is making
them reproducible and transferable, i.e. we should get the same results using the
pipeline on the same data, irrespective of which compute/system we analysed them
on. To achieve this, we will work with contained software environments, using
`conda` (please install [Anaconda](https://www.anaconda.com/products/individual)
if you haven't yet, so we can make use of this feature).

During the lectures, we will build our pipeline from scratch, starting with how
to write snakemake rules. In the end, we will have built an analysis pipeline for
RNAseq including alignment, quality control and differential expression analysis.
There is also a [homework](HOMEWORK.md) assignment that will ask you to extend
the pipeline, evaluate your results and learn how to generate an analysis report.

## Directory structure
The files you see in this directory are either part of the pipeline
([Snakemake](Snakemake),[scripts](scripts), [envs](envs)), or contain the example
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
not have one, give [Atom](https://atom.io/) a try!

**2. Install Snakemake**
Following the recommodations by the
[snakemake developers](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html),
we will first install _mamba_, a robuster and faster version of the default
conda package manager that is shipped with Anaconda.

- Open your terminal app (aka commad line, shell)
- Type the following command:

```bash
conda install -n base -c conda-forge mamba
```
NB: this might take a while. Initially you should see:
```bash
Collecting package metadata (current_repodata.json): done
Collecting package metadata (repodata.json): done
Solving environment: done
```
At some point it will ask you `Proceed ([y]/n)?`. Type `y` and hit enter.

If it successfully finishes, you will see this on your screen:
```bash
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
```

We will then use mamba to create a software environment specific for our
analysis. At the moment, all it needs to contain is snakemake itself. To create
this environment type:

```bash
mamba create -c conda-forge -c bioconda -n snakemake-rnaseq snakemake
```

As above, this might take a while and it will ask you again
`Proceed ([y]/n)?`. Type `y` and hit enter. Once it finishes (if successful) you
will see:

```bash
# To activate this environment, use
#
#     $ conda activate snakemake-rnaseq
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```

**3. Activate the snakemake environment**
- type the following to activate the environment
```bash
conda activate snakemake-rnaseq
```
- test the environment by typing
```bash
snakemake --version
```
You should see: `6.10.0`

**4. Move into the rnaseq-tutorial directory**
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

- type: `snakemake --use-conda --cores 1`
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
