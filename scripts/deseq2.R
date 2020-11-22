log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#################
## libraries ####
#################
library("DESeq2")
library("tidyverse")


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

dds <- readRDS(snakemake@input[["dds"]])

################
## analysis ####
################

## 1. Model fit ####
# Generate named coefficients need for apeglm lfcShrink
comparison <- snakemake@wildcards[["contrast"]]
coef <- paste("condition", comparison, sep="_")

# Relevel for reference to second element in contrasts
dds$condition <- relevel(dds$condition, "control")
dds <- nbinomWaldTest(dds)

## 2. Process results ####
# Extract coefficient specific results
res <- results(dds, name=coef, parallel=parallel)

# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, coef=coef, res=res, type='apeglm')

# add gene names to results object
res$gene_name <- mcols(dds)$symbol

# sort by p-value
res <- res[order(res$padj),]

## 3. Summarise results ####
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
