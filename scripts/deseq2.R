log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#################
## libraries ####
#################
library("DESeq2")
library("biomaRt")
library("tidyverse")
library("cowplot")

#################
## functions ####
#################
genes_de <- function(deset, thrP=0.05, thrLog2FC=log2(1.5),
                     direction=c('up', 'down', 'any')) {
    tmp <- deset %>%
        as.data.frame %>%
        rownames_to_column(var="gene_id")
    if (direction == "up") {
        tmp <- tmp %>%
            dplyr::filter(padj < thrP,
                          log2FoldChange > thrLog2FC)
    } else if (direction == "down") {
        tmp <- tmp %>%
            dplyr::filter(padj < thrP,
                          log2FoldChange < thrLog2FC)
    } else if (direction == "any") {
        tmp <- tmp %>%
            dplyr::filter(padj < thrP)
    } else {
        stop(direction, "is not a valid option for direction")
    }
    tmp
}

save_up_down <- function(res, snakemake) {
    up  <- genes_de(res, direction="up")
    down  <- genes_de(res, direction="down")
    write_csv(up, snakemake@output[['up']])
    write_csv(down, snakemake@output[['down']])
    return(list(up=up$gene_id, down=down$gene_id))
}

############
## data ####
############
message("Reading counts")
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene",
                  check.names=FALSE)

message("Reading sample file")
coldata <- read.table(snakemake@input[["samples"]], header=TRUE,
                      row.names="sample", check.names=FALSE, sep="\t")

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

# Annotate by gene names
if (snakemake@params[["organism"]] == "mouse") {
    ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
} else if(snakemake@params[["organism"]] == "human") {
    ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
} else {
    stop("No annotation (biomart) specified for organism:",
         snakemake@params[["organism"]])
}

genemap <- getBM(attributes=c('ensembl_gene_id', "external_gene_name"),
       filters = 'ensembl_gene_id',
       values = rownames(dds),
       mart = ensembl) %>%
    dplyr::rename(gene_id = ensembl_gene_id,
           symbol = external_gene_name) %>%
    distinct()

featureData <- tibble(gene_id=rownames(dds)) %>%
    left_join(genemap, by="gene_id") %>%
    mutate(symbol=case_when(is.na(symbol) ~ gene_id,
                            TRUE ~ symbol)) %>%
    select(symbol)
mcols(dds) <- DataFrame(mcols(dds), featureData)

## 2. Model fit ####
# Generate named coefficients need for apeglm lfcShrink
elements <- snakemake@params[["contrast"]]
comparison <- paste(elements[1], "vs", elements[2], sep="_")
coef <- paste("condition", elements[1], "vs", elements[2], sep="_")

# Relevel for reference to second element in contrasts
dds$condition <- relevel(dds$condition, elements[2])
dds <- nbinomWaldTest(dds)

## 3. Process results ####
# Extract coefficient specific results
res <- results(dds, name=coef, parallel=parallel)

# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, coef=coef, res=res, type='apeglm')

# add gene names to results object
res$gene_name <- mcols(dds)$symbol

# sort by p-value
res <- res[order(res$padj),]

## 4. Summarise results ####
## a) All genes for all groups ####
res_format <- res %>%
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    as_tibble() %>%
    rename_at(vars(-gene_id, -gene_name), ~ paste0(., "_", comparison))

# normalised expression values
rld <- rlog(dds, blind = FALSE)
deg_genes <- assay(rld)

combined <- deg_genes %>%
    as.data.frame() %>%
    rownames_to_column(var="gene_id") %>%
    as_tibble() %>%
    right_join(res_format, by="gene_id") %>%
    dplyr::select(gene_id, gene_name, everything())

write_delim(combined, snakemake@output[["table"]], delim="\t")

## b) Up/Down genes ####
genes_up_down <- save_up_down(res=res, snakemake=snakemake)

## 4. Visualise results ####
# ma plot
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

# pca plot
# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(counts, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition))
p <- p + geom_point() +
    scale_color_brewer(type="qual", palette="Dark2") +
    labs(x=paste0("PC1: ", percentVar[1], "% variance"),
         y=paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.justification = 0)
ggsave(plot=p, height=4.5, width=7.5, filename = snakemake@output[["pca_plot"]])
