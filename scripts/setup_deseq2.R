log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#################
## libraries ####
#################
library("DESeq2")
library("biomaRt")
library("tidyverse")

############
## data ####
############
message("Reading counts")
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene",
                  check.names=FALSE)

message("Reading sample file")
coldata <- read.table(snakemake@input[["samples"]], header=TRUE,
                      row.names="sample", check.names=FALSE, sep="\t")

message("Reading annotation file")
annotation <- read.table(snakemake@input[["annotation"]], header=TRUE,
                       check.names=FALSE, sep="\t")

message("Getting experimental design")
design <- as.formula(snakemake@params[["design"]])

# colData and countData must have the same sample order
if (nrow(coldata) != ncol(cts)) {
    stop("Number of samples in sample sheet and number of samples in counts",
         "matrix is not the same")
}
cts <- cts[,match(rownames(coldata),colnames(cts))]
if (any(c("control", "Control", "CONTROL") %in% levels(coldata$condition))) {
    if ("control" %in% levels(coldata$condition)) {
        coldata$condition <- relevel(coldata$condition, "control" )
    } else if ("Control" %in% levels(coldata$condition)) {
        coldata$condition <- relevel(coldata$condition, "Control" )
    } else {
        coldata$condition <- relevel(coldata$condition, "CONTROL" )
    }
}

################
## analysis ####
################

## 1. Generate and annotate DESeq object ####
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=design)

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 1,]

# normalization and preprocessing
dds <- DESeq(dds)

# Remove build number on ENS gene id
rownames(dds) <- gsub("\\.\\d*", "", rownames(dds))

genemap <- annotation %>%
    distinct()

featureData <- tibble(gene_id=rownames(dds)) %>%
    left_join(genemap, by="gene_id") %>%
    mutate(symbol=case_when(is.na(symbol) ~ gene_id,
                            TRUE ~ symbol)) %>%
    select(symbol)

mcols(dds) <- DataFrame(mcols(dds), featureData)

saveRDS(dds, snakemake@output[["dds"]])
