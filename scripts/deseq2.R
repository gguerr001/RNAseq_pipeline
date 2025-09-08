suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(pheatmap)
  library(openxlsx)
  library(dplyr)
  library(stringr)
  library(data.table)
})

option_list <- list(
  make_option(c("--counts"), type="character"),
  make_option(c("--samples"), type="character"),
  make_option(c("--contrasts"), type="character"),
  make_option(c("--outdir"), type="character", default="results/deseq2")
)
opt <- parse_args(OptionParser(option_list=option_list))

counts <- fread(opt$counts)
rownames_counts <- counts$Geneid
counts$Geneid <- NULL
mat <- as.matrix(counts)
rownames(mat) <- rownames_counts

samples <- fread(opt$samples)
rownames(samples) <- samples$sample
samples$sample <- NULL

stopifnot(all(colnames(mat) %in% rownames(samples)))
# reorder colData to match count matrix
samples <- samples[colnames(mat), , drop=FALSE]

# DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=mat,
                              colData=samples,
                              design=~ condition)
dds <- dds[rowSums(counts(dds)) > 1, ]

dds <- DESeq(dds)

# Normalized counts table
norm <- counts(dds, normalized=TRUE)
write.table(data.frame(gene=rownames(norm), norm),
            file=file.path(opt$outdir, "normalized_counts.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# PCA plot (VST)
vsd <- vst(dds, blind=TRUE)
plotdata <- as.data.frame(plotPCA(vsd, intgroup=c("condition"), returnData=TRUE))
percentVar <- round(100 * attr(plotdata, "percentVar"))
p <- ggplot(plotdata, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_minimal(base_size = 12)
ggsave(file.path(opt$outdir, "PCA.pdf"), p, width=6, height=5)

# Heatmap of top 50 most variable genes
rv <- rowVars(assay(vsd))
top <- order(rv, decreasing=TRUE)[seq_len(min(50, length(rv)))]
mat_scaled <- assay(vsd)[top, ]
mat_scaled <- t(scale(t(mat_scaled)))
ann <- data.frame(condition = colData(vsd)$condition)
rownames(ann) <- colnames(vsd)
pdf(file.path(opt$outdir, "heatmap_top50.pdf"), width=6, height=8)
pheatmap(mat_scaled, annotation_col = ann, show_rownames = FALSE)
dev.off()

# Differential expression for each contrast
contrasts_df <- fread(opt$contrasts)
for (i in seq_len(nrow(contrasts_df))){
  case <- contrasts_df$case[i]
  control <- contrasts_df$control[i]
  message(sprintf("Running contrast: %s vs %s", case, control))

  coef_name <- paste0("condition_", case, "_vs_", control)
  res <- NULL
  if (coef_name %in% resultsNames(dds)){
    # preferred: lfcShrink with coef if available
    res <- lfcShrink(dds, coef=coef_name, type="apeglm")
  } else {
    # fallback: standard results with contrast
    res <- results(dds, contrast=c("condition", case, control))
  }
  resOrdered <- as.data.frame(res[order(res$pvalue),])
  resOrdered$gene <- rownames(resOrdered)

  # merge normalized counts
  normdf <- as.data.frame(norm)
  normdf$gene <- rownames(normdf)
  final <- merge(resOrdered, normdf, by="gene", all.x=TRUE)

  # save TSV + XLSX
  outprefix <- file.path(opt$outdir, paste0(case, "_vs_", control))
  write.table(final, paste0(outprefix, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  openxlsx::write.xlsx(final, paste0(outprefix, ".xlsx"))
}