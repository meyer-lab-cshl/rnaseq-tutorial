import pandas as pd

print(snakemake.params.strand)
def get_column(strandedness):
    if pd.isnull(strandedness) or strandedness == "none":
        return 1 #non stranded protocol
    elif strandedness == "yes":
        return 2 #3rd column
    elif strandedness == "reverse":
        return 3 #4th column, usually for Illumina truseq
    else:
        raise ValueError(("'strandedness' column should be empty or have the "
                          "value 'none', 'yes' or 'reverse', instead has the "
                          "value {}").format(repr(strandedness)))

counts = [pd.read_table(f, index_col=0, usecols=[0, get_column(snakemake.params.strand)],
          header=None, skiprows=4) for f in snakemake.input]

for t, sample in zip(counts, snakemake.params.samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
if len(set(matrix.columns)) != len(matrix.columns):
    matrix = matrix.groupby(matrix.columns, axis=1).sum()
matrix.to_csv(snakemake.output[0], sep="\t")
