# Usage: Rscript deseq2_analysis.R <count_matrix> <metadata_csv> <output_dir>
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
})

# -----------------------------
# Args
# -----------------------------
args       <- commandArgs(trailingOnly = TRUE)
count_file <- args[1]
meta_file  <- args[2]
out_dir    <- args[3]

# -----------------------------
# Load data
# -----------------------------
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))  # DESeq2 requires integers

meta <- read.csv(meta_file, row.names = 1)

# Ensure sample order matches between counts and metadata
if (!all(colnames(counts) %in% rownames(meta))) {
    stop("Sample names in count matrix do not match metadata. Check your CSV files.")
}
meta <- meta[colnames(counts), , drop = FALSE]
meta$condition <- factor(meta$condition)

# -----------------------------
# DESeq2
# -----------------------------
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = meta,
    design    = ~condition
)

# Filter lowly expressed genes (at least 10 counts across all samples)
dds <- dds[rowSums(counts(dds)) >=0, ]

dds <- DESeq(dds)

# Get contrast: second level vs first level (treated vs control)
contrast <- c("condition",
              levels(meta$condition)[2],
              levels(meta$condition)[1])

res <- results(dds, contrast = contrast, alpha = 0.01)
res <- res[order(res$padj), ]

# Save results
res_df <- as.data.frame(res)
res_df <- cbind(gene = rownames(res_df), res_df)
write.csv(res_df, file.path(out_dir, "deseq2_results.csv"), row.names = FALSE)

cat(sprintf("Total genes tested: %d\n", nrow(res_df)))
cat(sprintf("Significant (padj < 0.01): %d\n", sum(res$padj < 0.01, na.rm = TRUE)))
cat(sprintf("Upregulated: %d\n",   sum(res$padj < 0.01 & res$log2FoldChange >= 1, na.rm = TRUE)))
cat(sprintf("Downregulated: %d\n", sum(res$padj < 0.01 & res$log2FoldChange <= -1, na.rm = TRUE)))

# -----------------------------
# PCA plot
# -----------------------------
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle("PCA - Variance Stabilized Counts") +
    theme_bw()

ggsave(file.path(out_dir, "pca_plot.pdf"), pca_plot, width = 7, height = 5)

# -----------------------------
# Volcano plot
# -----------------------------
vol_df <- as.data.frame(res)
vol_df$gene      <- rownames(vol_df)
vol_df$sig       <- ifelse(is.na(vol_df$padj), "NA",
                    ifelse(vol_df$padj < 0.01 & abs(vol_df$log2FoldChange) >= 1,
                           "Significant", "Not significant"))
vol_df$label     <- ifelse(vol_df$sig == "Significant", vol_df$gene, NA)

volcano_plot <- ggplot(vol_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20) +
    scale_color_manual(values = c(
        "Significant"     = "red",
        "Not significant" = "grey50",
        "NA"              = "grey80"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted p-value") +
    ggtitle(paste("Volcano Plot:", levels(meta$condition)[2],
                  "vs", levels(meta$condition)[1])) +
    theme_bw() +
    theme(legend.title = element_blank())

ggsave(file.path(out_dir, "volcano_plot.pdf"), volcano_plot, width = 8, height = 6)

cat("✅ All outputs saved.\n")
