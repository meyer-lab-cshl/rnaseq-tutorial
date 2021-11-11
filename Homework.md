# Homework

During the course, we developed a pipeline for RNAseq analysis, from generating
genome indices and quality control, to trimming, aligning and differential
expression analysis. This week's homework includes creating additional
rules for QC and visualising the results of the expression expression
analysis using PCA.

## Quality control
Look at the documentation of reseqc: http://rseqc.sourceforge.net/
Find the usage information on how to figure out the strandedness of your library.
Include a rule calling that command in your workflow and visualise via
your multiqc rule.

## PCA
Use the pca.R script in the scripts folder to visualise the expression data
in a PCA plot. As input it requires the 'deseq2/all.rds' output from the
setup_de rule named as dds. As output it requires the name 'pca_plot'; you
can choose the output file name, but it should end in 'pdf' (have a look at the
script file if you are interested what's happening under the hood.). The output
name in the rule should follow the way we have called the table output in the
'deseq2' rule. Use the deseq2.yaml as conda environment.

## Differentially expressed genes
Have a look at the output files of 'deseq2' rule. Use the command line to figure
out how many up and down-regulated genes there are (no need to include in
snakemake).

## Final report
Snakemake has the neat functionality of generating a final html report. Have
a look at this link: https://snakemake.readthedocs.io/en/stable/snakefiles/reporting.html
and generate a final html report of your pipeline (make sure you've got the
latest snakemake installed: version 6.10; if you got an earlier version,
there might be a bug in the snakemake report function, throwing this error:
`NameError: name 'contains_wildcard' is not defined`).

Send your final snakemake file, your html report and the number of
differentially expressed genes for grading.
