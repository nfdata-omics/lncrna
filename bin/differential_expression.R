#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(optparse)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
})

# Parse arguments
option_list <- list(
    make_option(c("--counts"), type = "character", help = "Count matrix file"),
    make_option(c("--gene_info"), type = "character", help = "Gene information file"),
    make_option(c("--design"), type = "character", help = "Design matrix file"),
    make_option(c("--method"), type = "character", default = "DESeq2", help = "DE method: DESeq2, edgeR, or limma [default: %default]"),
    make_option(c("--output"), type = "character", default = "DE", help = "Output prefix"),
    make_option(c("--padj_cutoff"), type = "double", default = 0.05, help = "Adjusted p-value cutoff"),
    make_option(c("--lfc_cutoff"), type = "double", default = 1, help = "Log2 fold change cutoff")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Method:", opt$method, "\n")

# Load DE libraries conditionally
if (opt$method == "DESeq2") {
    suppressPackageStartupMessages(library(DESeq2))
} else if (opt$method %in% c("edgeR", "limma")) {
    if (opt$method == "edgeR" && !requireNamespace("edgeR", quietly = TRUE)) {
        stop("Package 'edgeR' is required but not installed.")
    }
    if (opt$method == "limma" && !requireNamespace("limma", quietly = TRUE)) {
        stop("Package 'limma' is required but not installed.")
    }

    if (requireNamespace("edgeR", quietly = TRUE)) suppressPackageStartupMessages(library(edgeR))
    if (requireNamespace("limma", quietly = TRUE)) suppressPackageStartupMessages(library(limma))
}

# Read input data
cat("Reading inputs...\n")
counts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
gene_info <- read.table(opt$gene_info, header = TRUE, row.names = 1, sep = "\t", quote = "")
cat("Gene info loaded:", nrow(gene_info), "genes. Columns:", paste(colnames(gene_info), collapse=", "), "\n")

# Read design file (expecting CSV)
design <- read.csv(opt$design, header = TRUE, stringsAsFactors = TRUE)

# Standardize design column names
if ("sample" %in% colnames(design)) {
    colnames(design)[colnames(design) == "sample"] <- "SampleID"
} else if (!"SampleID" %in% colnames(design)) {
    # If SampleID is not present and sample is not present, assume first column is ID
    warning("Neither 'SampleID' nor 'sample' column found in design file. Using first column as sample ID.")
    colnames(design)[1] <- "SampleID"
}

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

# Helper function to safely add gene info
add_gene_info <- function(res_df, gene_info_df) {
    # Check if gene_info_df is valid
    if (is.null(gene_info_df) || nrow(gene_info_df) == 0) {
        warning("gene_info is empty or NULL. Skipping annotation.")
        res_df$gene_name <- rownames(res_df)
        res_df$gene_type <- "unknown"
        return(res_df)
    }

    # Use match for safer indexing
    idx <- match(rownames(res_df), rownames(gene_info_df))

    # Safely retrieve gene_name
    if ("gene_name" %in% colnames(gene_info_df)) {
        res_df$gene_name <- gene_info_df$gene_name[idx]
    } else {
        # Try to find a column that looks like gene name
        name_col <- grep("name|symbol", colnames(gene_info_df), ignore.case = TRUE, value = TRUE)
        if (length(name_col) > 0) {
            warning(paste("'gene_name' column not found. Using", name_col[1]))
            res_df$gene_name <- gene_info_df[[name_col[1]]][idx]
        } else {
            warning("'gene_name' column not found in gene_info. Using gene IDs.")
            res_df$gene_name <- rownames(res_df)
        }
    }

    # Safely retrieve gene_type
    if ("gene_type" %in% colnames(gene_info_df)) {
        res_df$gene_type <- gene_info_df$gene_type[idx]
    } else {
        # Try to find a column that looks like type or biotype
        type_col <- grep("type|biotype", colnames(gene_info_df), ignore.case = TRUE, value = TRUE)
        if (length(type_col) > 0) {
             warning(paste("'gene_type' column not found. Using", type_col[1]))
             res_df$gene_type <- gene_info_df[[type_col[1]]][idx]
        } else {
            warning("'gene_type' column not found in gene_info. Setting to 'unknown'.")
            res_df$gene_type <- "unknown"
        }
    }

    # Fill NAs (in case of no match)
    res_df$gene_name[is.na(res_df$gene_name)] <- rownames(res_df)[is.na(res_df$gene_name)]
    res_df$gene_type[is.na(res_df$gene_type)] <- "unknown"

    return(res_df)
}

# Run differential expression based on method
if (opt$method == "DESeq2") {
    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = design, design = ~condition)

    # Filter low count genes
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    cat("Genes after filtering:", nrow(dds), "\n")

    # Run DESeq2
    dds <- DESeq(dds)
    cat("DESeq2 analysis finished.\n")

    # Get results
    cat("Extracting results...\n")
    res_de <- results(dds, alpha = opt$padj_cutoff)
    res <- as.data.frame(res_de)
    cat("Results extracted.\n")

    # Add gene info
    cat("Adding gene info...\n")
    res <- add_gene_info(res, gene_info)
    cat("Gene info added.\n")

    # Normalized counts
    cat("Extracting normalized counts...\n")
    norm_counts <- counts(dds, normalized = TRUE)
    cat("Normalized counts extracted.\n")

    # Save DESeq2 object
    cat("Saving DESeq2 object...\n")
    saveRDS(dds, file = paste0(opt$output, ".dds.rds"))
    cat("DESeq2 object saved.\n")

    # VST transformation (faster and more robust than rlog for larger datasets or when rlog fails)
    # Use tryCatch to fall back or report error
    tryCatch({
        rld <- vst(dds, blind = FALSE)
        cat("VST transformation successful.\n")
    }, error = function(e) {
        warning("VST failed, trying rlog (slower): ", e$message)
        rld <<- rlog(dds, blind = FALSE)
    })

    # Dispersion plot
    tryCatch({
        # Extract dispersion data for ggplot
        # Based on DESeq2::plotDispEsts logic
        dds_mcols <- mcols(dds)
        df_disp <- data.frame(
            baseMean = dds_mcols$baseMean,
            dispGeneEst = dds_mcols$dispGeneEst,
            dispersion = dds_mcols$dispersion,
            dispFit = dds_mcols$dispFit
        )
        df_disp <- df_disp[df_disp$baseMean > 0 & !is.na(df_disp$dispersion), ]

        p <- ggplot(df_disp, aes(x = baseMean, y = dispersion)) +
            geom_point(aes(y = dispGeneEst), color = "black", alpha = 0.5, size = 0.5) +
            geom_point(color = "red", alpha = 0.5, size = 0.5) +
            scale_x_log10() +
            scale_y_log10() +
            labs(title = "Dispersion Estimates", x = "Mean of Normalized Counts", y = "Dispersion") +
            theme_bw()

        if (!all(is.na(df_disp$dispFit))) {
             p <- p + geom_line(aes(y = dispFit), color = "red")
        }

        png(paste0(opt$output, ".dispersion_plot.png"), width = 9, height = 7, units = "in", res = 150, type = "Xlib")
        print(p)
        dev.off()

        cat("Dispersion plot saved.\n")
    }, error = function(e) { warning("Failed to generate dispersion plot: ", e$message) })

    # PCA plot
    cat("Generating PCA plot...\n")
    tryCatch({
        pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
        percentVar <- round(100 * attr(pca_data, "percentVar"))

        p <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
            geom_point(size = 4) +
            xlab(paste0("PC1: ", percentVar[1], "% variance")) +
            ylab(paste0("PC2: ", percentVar[2], "% variance")) +
            ggtitle("PCA Plot - Sample Relationships") +
            theme_bw(base_size = 14) +
            theme(legend.position = "bottom")

        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(size = 3)
        } else {
            p <- p + geom_text(size = 3, vjust = 1.5)
        }

        png(paste0(opt$output, ".pca_plot.png"), width = 9, height = 7, units = "in", res = 150, type = "Xlib")
        print(p)
        dev.off()
        cat("PCA plot saved.\n")
    }, error = function(e) { warning("Failed to generate PCA plot: ", e$message) })

    # MA plot
    cat("Generating MA plot...\n")
    tryCatch({
        # Recreate MA plot using ggplot
        df_ma <- as.data.frame(res_de)
        df_ma <- df_ma[!is.na(df_ma$baseMean) & !is.na(df_ma$log2FoldChange), ]
        df_ma$significant <- ifelse(df_ma$padj < opt$padj_cutoff, "Yes", "No")
        df_ma$significant[is.na(df_ma$significant)] <- "No"

        p <- ggplot(df_ma, aes(x = baseMean, y = log2FoldChange, color = significant)) +
            geom_point(size = 0.8, alpha = 0.6) +
            scale_x_log10() +
            scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
            geom_hline(yintercept = 0, color = "firebrick", linetype = "dashed") +
            ylim(-5, 5) +
            labs(title = "MA Plot", x = "Mean of Normalized Counts", y = "Log2 Fold Change") +
            theme_bw() +
            theme(legend.position = "none")

        png(paste0(opt$output, ".ma_plot.png"), width = 9, height = 7, units = "in", res = 150, type = "Xlib")
        print(p)
        dev.off()
        cat("MA plot saved.\n")
    }, error = function(e) { warning("Failed to generate MA plot: ", e$message) })

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
    res <- add_gene_info(res, gene_info)

    # Normalized counts (CPM)
    norm_counts <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)

    # For heatmap, use log2(CPM+1)
    mat_for_heatmap <- cpm(y, normalized.lib.sizes = TRUE, log = TRUE)

    # MDS plot (similar to PCA)
    tryCatch({
        # Recreate MDS plot with ggplot
        # Calculate MDS
        mds <- limma::plotMDS(y, plot = FALSE)
        df_mds <- data.frame(
            Dim1 = mds$x,
            Dim2 = mds$y,
            SampleID = design$SampleID,
            Condition = design$condition
        )

        p <- ggplot(df_mds, aes(x = Dim1, y = Dim2, color = Condition, label = SampleID)) +
            geom_point(size = 4) +
            labs(title = "MDS Plot - Sample Relationships", x = "Dimension 1", y = "Dimension 2") +
            theme_bw(base_size = 14) +
            theme(legend.position = "bottom")

        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(size = 3)
        } else {
            p <- p + geom_text(size = 3, vjust = 1.5)
        }

        png(paste0(opt$output, ".pca_plot.png"), width = 9, height = 7, units = "in", res = 150, type = "Xlib")
        print(p)
        dev.off()
        cat("MDS plot saved.\n")
    }, error = function(e) { warning("Failed to generate MDS plot: ", e$message) })

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
    res <- add_gene_info(res, gene_info)

    # Normalized counts (CPM)
    norm_counts <- cpm(y, normalized.lib.sizes = TRUE, log = FALSE)

    # For heatmap
    mat_for_heatmap <- v$E

    # MDS plot
    tryCatch({
        # Recreate MDS plot with ggplot
        # Calculate MDS
        mds <- limma::plotMDS(y, plot = FALSE)
        df_mds <- data.frame(
            Dim1 = mds$x,
            Dim2 = mds$y,
            SampleID = design$SampleID,
            Condition = design$condition
        )

        p <- ggplot(df_mds, aes(x = Dim1, y = Dim2, color = Condition, label = SampleID)) +
            geom_point(size = 4) +
            labs(title = "MDS Plot - Sample Relationships", x = "Dimension 1", y = "Dimension 2") +
            theme_bw(base_size = 14) +
            theme(legend.position = "bottom")

        if (requireNamespace("ggrepel", quietly = TRUE)) {
            p <- p + ggrepel::geom_text_repel(size = 3)
        } else {
            p <- p + geom_text(size = 3, vjust = 1.5)
        }

        png(paste0(opt$output, ".pca_plot.png"), width = 9, height = 7, units = "in", res = 150, type = "Xlib")
        print(p)
        dev.off()
        cat("MDS plot saved.\n")
    }, error = function(e) { warning("Failed to generate MDS plot: ", e$message) })

} else {
    stop("Invalid DE method. Must be DESeq2, edgeR, or limma")
}

# Write results
res_ordered <- res[order(res$padj), ]
write.table(res_ordered, file = paste0(opt$output, ".de_results.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(norm_counts, file = paste0(opt$output, ".normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# Generate common plots

# Volcano plot
tryCatch({
    df_plot <- as.data.frame(res)
    df_plot$significant <- ifelse(
        abs(df_plot$log2FoldChange) > opt$lfc_cutoff & df_plot$padj < opt$padj_cutoff,
        ifelse(df_plot$log2FoldChange > 0, "Up", "Down"),
        "NS"
    )

    # Count significant genes
    n_up <- sum(df_plot$significant == "Up", na.rm = TRUE)
    n_down <- sum(df_plot$significant == "Down", na.rm = TRUE)

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

    png(paste0(opt$output, ".volcano_plot.png"), width = 10, height = 8, units = "in", res = 150, type = "Xlib")
    print(p)
    dev.off()
    cat("Volcano plot saved.\n")
}, error = function(e) { warning("Failed to generate Volcano plot: ", e$message) })

# Heatmap of top DE genes
tryCatch({
    top_genes <- head(rownames(res_ordered[!is.na(res_ordered$padj), ]), 50)

    if (length(top_genes) > 0) {
        mat <- mat_for_heatmap[top_genes, , drop = FALSE]

        # Remove zero variance rows to avoid scaling errors
        row_vars <- apply(mat, 1, var)
        mat <- mat[row_vars > 0, , drop = FALSE]

        if (nrow(mat) > 1) {
            # Scale rows for better visualization
            mat_scaled <- t(scale(t(mat)))

            # Annotation
            annotation_col <- data.frame(Condition = design$condition)
            rownames(annotation_col) <- design$SampleID

            png(paste0(opt$output, ".heatmap.png"), width = 10, height = 12, units = "in", res = 150, type = "Xlib")
            print(pheatmap(mat_scaled,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                annotation_col = annotation_col,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                main = "Top 50 DE Genes (Z-score)",
                fontsize_row = 8
            ))
            dev.off()
            cat("Heatmap saved.\n")
        } else {
            warning("Not enough genes with variance for heatmap.")
        }
    }
}, error = function(e) { warning("Failed to generate Heatmap: ", e$message) })
