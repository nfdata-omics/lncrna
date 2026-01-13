#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(DESeq2)
    library(edgeR)
    library(limma)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(ggrepel)
})

# Parse arguments
option_list <- list(
    make_option(c("--counts"), type = "character", help = "Count matrix file"),
    make_option(c("--gene_info"), type = "character", help = "Gene information file"),
    make_option(c("--design"), type = "character", help = "Design matrix file"),
    make_option(c("--method"), type = "character", default = "DESeq2", help = "DE method: DESeq2, edgeR, or limma [default: %default]"),
    make_option(c("--output"), type = "character", default = "DE", help = "Output prefix"),
    make_option(c("--cores"), type = "integer", default = 1, help = "Number of cores"),
    make_option(c("--padj_cutoff"), type = "double", default = 0.05, help = "Adjusted p-value cutoff"),
    make_option(c("--lfc_cutoff"), type = "double", default = 1, help = "Log2 fold change cutoff")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Method:", opt$method, "\n")

# Read input data
counts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
design <- read.table(opt$design, header = TRUE, sep = "\t", stringsAsFactors = TRUE)
gene_info <- read.table(opt$gene_info, header = TRUE, row.names = 1, sep = "\t", quote = "")

# Validate inputs
if (!all(design$SampleID %in% colnames(counts))) {
    stop("Not all samples in design file are present in count matrix")
}

# Ensure sample order matches
counts <- counts[, match(design$SampleID, colnames(counts))]

# Check if we have at least 2 conditions
if (length(unique(design$condition)) < 2) {
    stop("Need at least 2 conditions for differential expression analysis")
}

# Convert counts to integer matrix
counts <- round(as.matrix(counts))

cat("Samples:", ncol(counts), "\n")
cat("Genes:", nrow(counts), "\n")
cat("Conditions:", paste(unique(design$condition), collapse = ", "), "\n")

# Run differential expression based on method
if (opt$method == "DESeq2") {
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~condition)

    # Filter low count genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat("Genes after filtering:", nrow(dds), "\n")

    # Run DESeq2
    dds <- DESeq(dds, parallel = opt$cores > 1, BPPARAM = MulticoreParam(opt$cores))


    # Get results
    res <- results(dds, alpha = opt$padj_cutoff)
    res <- as.data.frame(res)

    # Add gene info
    res$gene_name <- gene_info[rownames(res), "gene_name"]
    res$gene_type <- gene_info[rownames(res), "gene_type"]

    # Normalized counts
    norm_counts <- counts(dds, normalized = TRUE)

    # Save DESeq2 object
    saveRDS(dds, file = paste0(opt$output, ".dds.rds"))

    # Regularized log transformation for plots
    rld <- rlog(dds, blind = FALSE)

    # Dispersion plot
    png(paste0(opt$output, ".dispersion_plot.png"), width = 800, height = 600, res = 120)
    plotDispEsts(dds, main = "Dispersion Estimates")
    dev.off()

    # PCA plot
    png(paste0(opt$output, ".pca_plot.png"), width = 900, height = 700, res = 120)
    pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pca_data, "percentVar"))
    print(ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) +
        geom_text_repel(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        ggtitle("PCA Plot - Sample Relationships") +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom"))
    dev.off()

    # MA plot
    png(paste0(opt$output, ".ma_plot.png"), width = 900, height = 700, res = 120)
    plotMA(res, ylim = c(-5, 5), main = "MA Plot - Differential Expression", alpha = opt$padj_cutoff)
    dev.off()

    mat_for_heatmap <- assay(rld)
} else if (opt$method == "edgeR") {
    # Create DGEList
    y <- DGEList(counts = counts, group = design$condition)

    # Filter low count genes
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat("Genes after filtering:", nrow(y), "\n")

    # Normalize
    y <- calcNormFactors(y)

    # Design matrix
    design_matrix <- model.matrix(~condition, data = design)

    # Estimate dispersion
    y <- estimateDisp(y, design_matrix)

    # Fit model
    fit <- glmQLFit(y, design_matrix)
    qlf <- glmQLFTest(fit, coef = 2)

    # Get results
    res <- topTags(qlf, n = Inf, adjust.method = "BH")$table
    res$padj <- res$FDR
    res$log2FoldChange <- res$logFC

    # Add gene info
    res$gene_name <- gene_info[rownames(res), "gene_name"]
    res$gene_type <- gene_info[rownames(res), "gene_type"]

    # Normalized counts (CPM)
    norm_counts <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)

    # For heatmap, use log2(CPM+1)
    mat_for_heatmap <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE)

    # MDS plot (similar to PCA)
    png(paste0(opt$output, ".pca_plot.png"), width = 900, height = 700, res = 120)
    plotMDS(y, labels = design$SampleID, col = as.numeric(design$condition), main = "MDS Plot - Sample Relationships")
    legend("topright", legend = levels(design$condition), col = 1:length(levels(design$condition)), pch = 16)
    dev.off()
} else if (opt$method == "limma") {
    # Create DGEList
    y <- DGEList(counts = counts, group = design$condition)

    # Filter low count genes
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    cat("Genes after filtering:", nrow(y), "\n")

    # Normalize
    y <- calcNormFactors(y)

    # Design matrix
    design_matrix <- model.matrix(~condition, data = design)

    # Voom transformation
    v <- voom(y, design_matrix, plot = FALSE)

    # Fit model
    fit <- lmFit(v, design_matrix)
    fit <- eBayes(fit)

    # Get results
    res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    res$padj <- res$adj.P.Val
    res$log2FoldChange <- res$logFC

    # Add gene info
    res$gene_name <- gene_info[rownames(res), "gene_name"]
    res$gene_type <- gene_info[rownames(res), "gene_type"]

    # Normalized counts (CPM)
    norm_counts <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)

    # For heatmap
    mat_for_heatmap <- v$E

    # MDS plot
    png(paste0(opt$output, ".pca_plot.png"), width = 900, height = 700, res = 120)
    plotMDS(y, labels = design$SampleID, col = as.numeric(design$condition), main = "MDS Plot - Sample Relationships")
    legend("topright", legend = levels(design$condition), col = 1:length(levels(design$condition)), pch = 16)
    dev.off()
} else {
    stop("Invalid DE method. Must be DESeq2, edgeR, or limma")
}

# Write results
res_ordered <- res[order(res$padj), ]
write.table(res_ordered, file = paste0(opt$output, ".de_results.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(norm_counts, file = paste0(opt$output, ".normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Generate common plots

# Volcano plot
df_plot <- as.data.frame(res)
df_plot$significant <- ifelse(
    abs(df_plot$log2FoldChange) > opt$lfc_cutoff & df_plot$padj < opt$padj_cutoff,
    ifelse(df_plot$log2FoldChange > 0, "Up", "Down"),
    "NS"
)

# Count significant genes
n_up <- sum(df_plot$significant == "Up", na.rm = TRUE)
n_down <- sum(df_plot$significant == "Down", na.rm = TRUE)

png(paste0(opt$output, ".volcano_plot.png"), width = 1000, height = 800, res = 120)
p <- ggplot(df_plot, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "blue", "NS" = "grey", "Up" = "red"), labels = c(paste("Down:", n_down), "NS", paste("Up:", n_up))) +
    geom_vline(xintercept = c(-opt$lfc_cutoff, opt$lfc_cutoff), linetype = "dashed", color = "grey30") +
    geom_hline(yintercept = -log10(opt$padj_cutoff), linetype = "dashed", color = "grey30") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom") +
    labs(
        title = "Volcano Plot - Differential Expression",
        subtitle = paste0("Cutoffs: |log2FC| > ", opt$lfc_cutoff, ", padj < ", opt$padj_cutoff),
        x = "log2 Fold Change",
        y = "-log10 adjusted p-value",
        color = "Regulation"
    )
print(p)
dev.off()

# Heatmap of top DE genes
top_genes <- head(rownames(res_ordered[!is.na(res_ordered$padj), ]), 50)

if (length(top_genes) > 0) {
    mat <- mat_for_heatmap[top_genes, , drop = FALSE]

    # Scale rows for better visualization
    mat_scaled <- t(scale(t(mat)))

    # Annotation
    annotation_col <- data.frame(Condition = design$condition)
    rownames(annotation_col) <- design$SampleID

    png(paste0(opt$output, ".heatmap.png"), width = 1000, height = 1200, res = 120)
    pheatmap(mat_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        annotation_col = annotation_col,
        color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
        main = "Top 50 Differentially Expressed Genes",
        fontsize_row = 8
    )
    dev.off()
}
