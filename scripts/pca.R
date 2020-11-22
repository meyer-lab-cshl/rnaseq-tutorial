log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

#################
## libraries ####
#################
library("DESeq2")
library("tidyverse")
library("cowplot")

############
## data ####
############

dds <- readRDS(snakemake@input[["dds"]])
if (!is.null(snakemake@params[['colorby']])) {
  color <- sym(snakemake@params[['colorby']])
} else {
  color <- sym('condition')
}
################
## analysis ####
################

# pca plot
# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(counts, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = !!color))
p <- p + geom_point() +
    scale_color_brewer(type="qual", palette="Dark2") +
    labs(x=paste0("PC1: ", percentVar[1], "% variance"),
         y=paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_cowplot() +
    theme(legend.position = "bottom",
          legend.justification = 0)
ggsave(plot=p, height=4.5, width=7.5, filename = snakemake@output[["pca_plot"]])
