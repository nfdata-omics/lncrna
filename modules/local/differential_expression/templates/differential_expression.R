#!/usr/bin/env Rscript

# DGE edgeR script for lncRNAs
# Based on: DGE edgeR script for miRNAs (Karla Ruiz - IU2 OMICS HT, 2026-02-23)

#----- Load libraries ----
suppressMessages(library(tidyverse))
suppressMessages(library(tibble))
suppressMessages(library(edgeR))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(RColorBrewer))
suppressMessages(library(gplots))
suppressMessages(library(openxlsx))

# ---- Nextflow template variables ----
lnc_counts_file <- "${count_matrix}"
metadata_file   <- "${design_file}"
contrast        <- "${contrast}"

# ---- Parameters (from task.ext / params) with defaults ----
LFC_cutoff        <- as.numeric("${LFC}"         ); if (is.na(LFC_cutoff))  LFC_cutoff  <- 0.5
FDR_cutoff        <- as.numeric("${FDR}"         ); if (is.na(FDR_cutoff))  FDR_cutoff  <- 0.05
min_count         <- as.integer("${min_count}"   ); if (is.na(min_count))   min_count   <- 3
min_total         <- as.integer("${min_total}"   ); if (is.na(min_total))   min_total   <- 15
correction_mode   <- "${correction_mode}"; if (correction_mode == "" || is.na(correction_mode)) correction_mode <- "none"
batch             <- "${batch}"
ruv_k             <- as.integer("${ruv_k}"       ); if (is.na(ruv_k))       ruv_k       <- 1
ruv_controls      <- "${ruv_controls}"; if (ruv_controls == "" || is.na(ruv_controls)) ruv_controls <- "auto"
top_genes         <- as.integer("${top_genes}"   ); if (is.na(top_genes))   top_genes   <- 100

contrast_lower <- tolower(contrast)

cat("lncRNA counts file :", lnc_counts_file, "\n")
cat("Metadata file      :", metadata_file,   "\n")
cat("Contrast           :", contrast,         "\n")
cat("LFC cutoff         :", LFC_cutoff,          "\n")
cat("FDR cutoff         :", FDR_cutoff,          "\n")
cat("min.count          :", min_count,    "\n")
cat("min.total.count    :", min_total,    "\n")

# ============================================================
#  FUNCTIONS
# ============================================================

# ---- 1. Import data ----
import_data <- function(counts_file, metadata_file) {

  raw <- read.table(counts_file,
                    header        = TRUE,
                    sep           = "\t",
                    row.names     = 1,
                    check.names   = FALSE,
                    quote         = "",
                    comment.char  = "")

  raw <- as.data.frame(raw)
  # Round to integers and validate (edgeR requires integer counts)
  raw <- round(raw)
  if (any(raw < 0, na.rm = TRUE)) {
    stop("Count matrix contains negative values.")
  }
  raw <- as.matrix(raw)
  storage.mode(raw) <- "integer"

  # Clean rownames (replace dots with dashes)
  rownames(raw) <- gsub("\\\\.", "-", rownames(raw))

  # METADATA
  metadata <- read.delim(metadata_file, sep = ",", row.names = 1)

  cat("Metadata loaded:", nrow(metadata), "samples\n")
  cat("Samples in metadata:", paste(rownames(metadata), collapse = ", "), "\n")

  # Validate: all metadata samples must be present in count matrix
  missing <- setdiff(rownames(metadata), colnames(raw))
  if (length(missing) > 0) {
    stop(paste("Samples in metadata not found in count matrix:", paste(missing, collapse = ", ")))
  }

  # Reorder count matrix columns to match metadata row order
  raw <- raw[, rownames(metadata)]

  return(list(rawCountTable = raw, metadata = metadata))
}


# ---- 2. Identify condition column in metadata ----
find_condition_column <- function(metadata, contrast_str) {
  conditions <- unlist(strsplit(tolower(contrast_str), "[:/]"))
  col_name <- NULL
  for (cond in conditions) {
    found <- sapply(metadata, function(col) any(grepl(cond, tolower(as.character(col)))))
    if (any(found)) {
      col_name <- names(metadata)[which(found)[1]]
      break
    }
  }
  if (is.null(col_name)) {
    stop("Could not find a metadata column matching the contrast levels: ", contrast_str)
  }
  return(col_name)
}


# ---- 3. Main DGE function (QLF pipeline) ----
perform_dge_lncrna <- function(rawCountTable, metadata, group, samples,
                               contrast_str,
                               LFC.cutoff, padj.cutoff,
                               min_count, min_total,
                               correction_mode = "none", batch_var = NULL,
                               ruv_k           = 1,
                               ruv_controls    = "auto",
                               top_n) {

  dge <- DGEList(counts = rawCountTable, group = group, samples = samples)

  # ----------------------------------------------------------
  # FILTERING
  # Low-count lncRNAs are filtered with relaxed thresholds
  # because lncRNAs have naturally lower expression than mRNAs.
  # filterByExpr uses min.count and min.total.count together.
  # ----------------------------------------------------------
  cat("Filtering low-expressed lncRNAs (min.count =", min_count,
      ", min.total.count =", min_total, ")...\n")
  keep <- filterByExpr(dge,
                       min.count       = min_count,
                       min.total.count = min_total)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  cat("lncRNAs retained after filtering:", nrow(dge), "\n")

  # ----------------------------------------------------------
  # NORMALIZATION (TMM)
  # ----------------------------------------------------------
  dge <- calcNormFactors(dge, method = "TMM")

  # Design matrix with explicit group levels (needed for makeContrasts)
  #design_mat <- model.matrix(~0 + group)
  #colnames(design_mat) <- levels(group)

  # ----------------------------------------------------------
  # BATCH CORRECTION — conditional on correction_mode
  # none  : no correction, simple ~0 + group model
  # batch : known batch variable added to model as covariate
  #         requires batch and condition to NOT be confounded
  # ruv   : RUVSeq correction for unknown or confounded batch
  #         estimates unwanted variation factors from data
  # ----------------------------------------------------------
  if (correction_mode == "batch") {
    if (is.null(batch_var)) {
      stop("--correction_mode batch requires --batch <column_name> in metadata.")
    }
    if (length(levels(batch_var)) < 2) {
      warning("Batch variable has only one level. Running without batch correction.")
      design_mat <- model.matrix(~0 + group)
      colnames(design_mat) <- levels(group)
    } else {
      design_mat <- model.matrix(~0 + group + batch_var)
      colnames(design_mat) <- c(levels(group),
                                paste0("batch_", levels(batch_var)[-1]))
      cat("Batch variable included in model. Levels:",
          paste(levels(batch_var), collapse = ", "), "\n")
    }

  } else if (correction_mode == "ruv") {
    suppressMessages(library(RUVSeq))
    cat("Running RUVSeq correction (k =", ruv_k, ")...\n")

    # Select negative control genes
    logcpm_pre <- cpm(dge, log = TRUE)

    if (ruv_controls == "auto") {
      cv_genes   <- apply(logcpm_pre, 1, function(x) sd(x) / (mean(abs(x)) + 0.1))
      n_ctrl     <- max(50, round(nrow(dge) * 0.1))
      ctrl_genes <- names(sort(cv_genes))[1:min(n_ctrl, length(cv_genes))]
      cat("RUV auto control genes:", length(ctrl_genes), "\n")
    } else {
      if (!file.exists(ruv_controls)) {
        stop(paste("RUV controls file not found:", ruv_controls))
      }
      ctrl_genes <- readLines(ruv_controls)
      ctrl_genes <- ctrl_genes[ctrl_genes != "" & ctrl_genes %in% rownames(dge)]
      if (length(ctrl_genes) == 0) {
        stop("None of the provided RUV control genes found in count matrix.")
      }
      cat("RUV user-provided control genes found:", length(ctrl_genes), "\n")
    }

    # Run RUVg
    count_mat <- as.matrix(dge\$counts)
    set_ruv    <- RUVg(count_mat, ctrl_genes, k = ruv_k)
    W_factors  <- set_ruv\$W
    cat("RUV factors estimated:", paste(colnames(W_factors), collapse = ", "), "\n")

    # Replace counts with RUV-normalized counts
    dge\$counts <- set_ruv\$normalizedCounts

    # Build design matrix including W factors
    W_df       <- as.data.frame(W_factors)
    design_mat <- model.matrix(~0 + group + ., data = W_df)
    colnames(design_mat) <- c(levels(group),
                              paste0("W_", seq_len(ruv_k)))
    cat("Design matrix includes RUV W factors.\n")

  } else {
    # none — default
    design_mat <- model.matrix(~0 + group)
    colnames(design_mat) <- levels(group)
    cat("No batch correction applied.\n")
  }

  # ----------------------------------------------------------
  # DISPERSION ESTIMATION — modern single-step approach
  # estimateDisp() replaces the three-step classic pipeline
  # (estimateGLMCommonDisp / TrendedDisp / TagwiseDisp).
  # It fits common + trended + tagwise dispersions internally
  # and is the current edgeR recommendation.
  # ----------------------------------------------------------
  dge <- estimateDisp(dge, design_mat, robust = TRUE)
  # robust = TRUE uses a weighted likelihood to protect against
  # outlier genes that could inflate common dispersion estimates.

  # ----------------------------------------------------------
  # GLOBAL QC PLOTS ( MDS + PCA + Heatmap of top variable genes)
  # ----------------------------------------------------------
  if (ncol(dge\$counts) > 2) {

    # MDS plot
    tryCatch({
      col_groups <- as.numeric(dge\$samples\$group)

      png("MDS.png", height = 7, width = 7, res = 600, units = "in")
      plotMDS(dge, col = col_groups, cex = 0.8, main = "MDS - lncRNA expression")
      legend("bottomright", legend = levels(dge\$samples\$group),
             col = 1:nlevels(dge\$samples\$group), pch = 20, horiz = TRUE, cex = 0.6)
      dev.off()

      pdf("MDS.pdf", height = 7, width = 7)
      plotMDS(dge, col = col_groups, cex = 0.8, main = "MDS - lncRNA expression")
      legend("bottomright", legend = levels(dge\$samples\$group),
             col = 1:nlevels(dge\$samples\$group), pch = 20, horiz = TRUE, cex = 0.6)
      dev.off()
      cat("MDS plot saved.\n")
    }, error = function(e) warning("MDS plot failed: ", e\$message))

    # Top variable lncRNAs
    tryCatch({
      logcounts  <- cpm(dge, log = TRUE)
      cpm_counts <- as.data.frame(cpm(dge))
      var_genes  <- apply(logcounts, 1, var)
      top_var_n  <- min(50, nrow(logcounts))
      select_var <- names(sort(var_genes, decreasing = TRUE))[1:top_var_n]
      hv_lcpm    <- logcounts[select_var, ]

      # PCA plot — computed from log-CPM of top 500 most variable lncRNAs (or all genes if fewer than 500 pass filtering)
      ntop    <- min(500, nrow(logcounts))
      top500  <- names(sort(var_genes, decreasing = TRUE))[1:ntop]
      pca_mat <- t(logcounts[top500, ])           # samples x genes
      pca_res <- prcomp(pca_mat, scale. = FALSE)  # log-CPM already on comparable scale

      pct_var <- round(100 * pca_res\$sdev^2 / sum(pca_res\$sdev^2), 1)

      df_pca <- data.frame(
        PC1       = pca_res\$x[, 1],
        PC2       = pca_res\$x[, 2],
        SampleID  = rownames(pca_res\$x),
        Condition = group
      )

      p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = Condition, label = SampleID)) +
        geom_point(size = 4) +
        geom_text(size = 3, vjust = -0.8, hjust = 0.5) +
        labs(
          title    = "PCA - lncRNA expression (log-CPM)",
          subtitle = paste0("Top ", ntop, " most variable lncRNAs"),
          x        = paste0("PC1: ", pct_var[1], "% variance"),
          y        = paste0("PC2: ", pct_var[2], "% variance")
        ) +
        theme_bw(base_size = 14) +
        theme(legend.position = "bottom")

      ggsave("PCA.png", p_pca, height = 7, width = 8, dpi = 600)
      ggsave("PCA.pdf", p_pca, height = 7, width = 8)

      ## HEATMAP

      mypalette <- brewer.pal(11, "RdYlBu")
      morecols  <- colorRampPalette(mypalette)
      col.cell  <- RColorBrewer::brewer.pal(8, "Set1")[as.numeric(group)]

      png("heatmap_global.png", height = 7, width = 7, res = 600, units = "in")
      heatmap.2(hv_lcpm, col = rev(morecols(50)), trace = "none",
                main = paste("Top", top_var_n, "variable lncRNAs (all samples)"),
                ColSideColors = col.cell, scale = "row",
                margins = c(10, 10), cexCol = 0.8, cexRow = 0.7)
      dev.off()

      pdf("heatmap_global.pdf", height = 7, width = 7)
      heatmap.2(hv_lcpm, col = rev(morecols(50)), trace = "none",
                main = paste("Top", top_var_n, "variable lncRNAs (all samples)"),
                ColSideColors = col.cell, scale = "row",
                margins = c(10, 10), cexCol = 0.8, cexRow = 0.7)
      dev.off()
    }, error = function(e) warning("Global heatmap failed: ", e\$message))

  } else {
    warning("Less than 3 samples — skipping global QC plots.")
    logcounts  <- cpm(dge, log = TRUE)
    cpm_counts <- as.data.frame(cpm(dge))
    var_genes  <- apply(logcounts, 1, var)
  }

  # ----------------------------------------------------------
  # BUILD CONTRASTS
  # ----------------------------------------------------------
  parts    <- strsplit(contrast_str, ":")[[1]]
  tests    <- unlist(strsplit(parts[1], "/"))
  controls <- unlist(strsplit(parts[2], "/"))

  contrast_strings <- sapply(seq_along(tests), function(i) {
    ctrl <- if (i <= length(controls)) controls[i] else controls[length(controls)]
    paste0(tests[i], "-", ctrl)
  })
  cat("Contrasts to test:", paste(contrast_strings, collapse = "; "), "\n")

  my.contrasts <- makeContrasts(contrasts = contrast_strings, levels = design_mat)
  colnames(my.contrasts) <- gsub("-", "_vs_", colnames(my.contrasts))

  # ----------------------------------------------------------
  # FIT MODEL — QLF pipeline
  # glmQLFit provides better control of type-I error vs glmFit/LRT
  # ----------------------------------------------------------
  cat("Fitting quasi-likelihood model (glmQLFit)...\n")
  fit_ql <- glmQLFit(dge, design_mat, robust = TRUE)
  # robust = TRUE uses a trimmed mean of the quasi-dispersion estimates

  # Also fit classic GLM for secondary LRT output (kept as reference)
  fit_lm <- glmFit(dge, design_mat)

  # ----------------------------------------------------------
  # PER-CONTRAST ANALYSIS
  # ----------------------------------------------------------
  results_list <- list()

  for (i in seq_len(ncol(my.contrasts))) {
    contrast_name <- colnames(my.contrasts)[i]
    contr         <- my.contrasts[, i]
    cat("\n--- Processing contrast:", contrast_name, "---\n")

    ctr_dir <- paste0("./DE_", contrast_name)
    dir.create(file.path(ctr_dir, "plots", "png"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(ctr_dir, "plots", "pdf"), recursive = TRUE, showWarnings = FALSE)
    png_dir <- file.path(ctr_dir, "plots", "png")
    pdf_dir <- file.path(ctr_dir, "plots", "pdf")

    # ---- PRIMARY: QLF test ----
    qlf <- glmQLFTest(fit_ql, contrast = contr)
    res_qlf <- topTags(qlf, n = nrow(dge), adjust.method = "BH")

    # ---- SECONDARY: LRT (for comparison / reference) ----
    lrt <- glmLRT(fit_lm, contrast = contr)
    res_lrt <- topTags(lrt, n = nrow(dge), adjust.method = "BH")

    # Summary tables
    sum_qlf <- summary(decideTests(qlf))
    sum_lrt <- summary(decideTests(lrt))

    res_DGE_lab <- as.data.frame(res_qlf)
    res_DGE_lab\$diffexpressed <- "NO"
    res_DGE_lab\$diffexpressed[res_DGE_lab\$logFC >=  LFC.cutoff & res_DGE_lab\$FDR <= padj.cutoff] <- "UP"
    res_DGE_lab\$diffexpressed[res_DGE_lab\$logFC <= -LFC.cutoff & res_DGE_lab\$FDR <= padj.cutoff] <- "DOWN"

    res_DGE_lab <- res_DGE_lab %>%
      mutate(lncRNA = rownames(res_DGE_lab)) %>%
      select(lncRNA, everything()) %>%
      as.data.frame() %>%
      { rownames(.) <- NULL; . }

    # Top N for heatmap
    if (nrow(res_DGE_lab) < 100) {
      topgenes <- nrow(res_DGE_lab)
    } else {
      topgenes <- 100
    }
    print(paste("Top variable genes in", contrast_name, ":", topgenes))

    # ---- Export results ----
    # QLF
    write.xlsx(res_DGE_lab,
               file = file.path(ctr_dir, paste0("DGE_QLF_", contrast_name, ".xlsx")),
               rowNames = FALSE)

    # LRT
    res_lrt_df <- as.data.frame(res_lrt)
    write.table(res_lrt_df,
                file = file.path(ctr_dir, paste0("DGE_LRT_", contrast_name, ".txt")),
                row.names = TRUE, quote = FALSE, sep = "\t")

    # Normalized CPM
    cpm_out <- as.data.frame(cpm(dge, normalized.lib.sizes = TRUE, log = FALSE))
    write.csv(cpm_out,
              file = file.path(ctr_dir, paste0("normalized_CPM_", contrast_name, ".csv")), row.names = TRUE, quote = FALSE)

    # Summary text
    out_summary <- file.path(ctr_dir, paste0("DGE_summary_", contrast_name, ".txt"))
    conn <- file(out_summary, open = "wt")
    cat(paste0(
      "Differential Expression Analysis — lncRNAs\n",
      "Contrast  : ", contrast_name, "\n",
      "Method    : edgeR QLF (primary) + LRT (secondary)\n",
      "LFC cutoff: ", LFC.cutoff, "\n",
      "FDR cutoff: ", padj.cutoff, "\n\n",
      "=== QLF summary ===\n"
    ), file = conn)
    for (j in seq_len(nrow(sum_qlf))) {
      cat(rownames(sum_qlf)[j], "\t", sum_qlf[j, ], "\n", file = conn)
    }
    cat("\n=== LRT summary (reference) ===\n", file = conn)
    for (j in seq_len(nrow(sum_lrt))) {
      cat(rownames(sum_lrt)[j], "\t", sum_lrt[j, ], "\n", file = conn)
    }
    close(conn)

    # ---- PLOTS ----
    plot_title <- tools::toTitleCase(gsub("_", " ", contrast_name))

    # Volcano plot (EnhancedVolcano)
    tryCatch({
      vplot <- EnhancedVolcano(res_qlf\$table,
                               lab      = rownames(res_qlf),
                               x        = "logFC",
                               y        = "FDR",
                               pCutoffCol  = "FDR",
                               labSize     = 3.0,
                               pCutoff     = padj.cutoff,
                               FCcutoff    = LFC.cutoff,
                               subtitle    = "edgeR QLF",
                               title       = plot_title,
                               legendPosition  = "right",
                               legendLabSize   = 10,
                               legendIconSize  = 3)
      ggsave(file.path(png_dir, paste0("volcanoplot_", contrast_name, ".png")),
             vplot, height = 7, width = 10, dpi = 600)
      ggsave(file.path(pdf_dir, paste0("volcanoplot_", contrast_name, ".pdf")),
             vplot, height = 7, width = 10)
    }, error = function(e) warning("Volcano plot failed: ", e\$message))

    # Heatmap for this contrast
    tryCatch({
      contr_levels  <- unlist(strsplit(contrast_name, "_vs_"))
      samples_keep  <- dge\$samples %>% filter(group %in% contr_levels)
      select_var_lnc <- names(sort(var_genes, decreasing = TRUE))[1:topgenes]
      hv_lcpm_ctr    <- logcounts[select_var_lnc,
                                  colnames(logcounts) %in% rownames(samples_keep),
                                  drop = FALSE]

      mypalette <- brewer.pal(11, "RdYlBu")
      morecols  <- colorRampPalette(mypalette)
      col.cell  <- RColorBrewer::brewer.pal(8, "Set1")[samples_keep\$group]

      png(file.path(png_dir, paste0("heatmap_", contrast_name, ".png")),
          height = 7, width = 7, res = 600, units = "in")
      heatmap.2(hv_lcpm_ctr, col = rev(morecols(50)), trace = "none",
                main     = paste("Top", topgenes, "variable lncRNAs\n", plot_title),
                ColSideColors = col.cell, scale = "row",
                margins  = c(10, 10), cexCol = 1, cexRow = 0.7)
      dev.off()

      pdf(file.path(pdf_dir, paste0("heatmap_", contrast_name, ".pdf")),
          height = 7, width = 7)
      heatmap.2(hv_lcpm_ctr, col = rev(morecols(50)), trace = "none",
                main     = paste("Top", topgenes, "variable lncRNAs\n", plot_title),
                ColSideColors = col.cell, scale = "row",
                margins  = c(10, 10), cexCol = 1, cexRow = 0.7)
      dev.off()
    }, error = function(e) warning("Contrast heatmap failed: ", e\$message))

    # Store all contrast results
    results_list[[contrast_name]] <- list(
      res_qlf       = res_qlf,
      res_lrt       = res_lrt,
      res_DGE_lab   = res_DGE_lab,
      sum_qlf       = sum_qlf,
      sum_lrt       = sum_lrt,
      dge           = dge,
      fit_ql        = fit_ql
    )
  }

  return(results_list)
}


# ============================================================
#  MAIN EXECUTION
# ============================================================
data_lnc  <- import_data(lnc_counts_file, metadata_file)
lnccounts <- data_lnc\$rawCountTable
metadata  <- data_lnc\$metadata

# Find condition column
contrast_lower <- tolower(contrast)
col_name <- find_condition_column(metadata, contrast_lower)
cat("Condition column identified:", col_name, "\n")

group   <- factor(tolower(metadata[[col_name]]))
levels(group) <- gsub(" ", "_", levels(group))
samples <- factor(rownames(metadata))

cat("Groups:", paste(levels(group), collapse = ", "), "\n")
cat("Samples:", nrow(metadata), "\n")

# Batch variable — solo si correction_mode es batch
batch_var <- NULL
if (correction_mode == "batch") {
  if (is.null(batch)) {
    stop("--correction_mode batch requires --batch <column_name>.")
  }
  if (!batch %in% colnames(metadata)) {
    stop(paste("Batch column '", batch, "' not found in metadata."))
  }
  batch_var <- factor(metadata[[batch]])
  print(paste("Batch variable:", batch,
              "| Levels:", paste(levels(batch_var), collapse = ", ")))
}

# Run analysis
dge_res <- perform_dge_lncrna(
  rawCountTable = lnccounts,
  metadata      = metadata,
  group         = group,
  samples       = samples,
  contrast_str  = contrast_lower,
  LFC.cutoff    = LFC_cutoff,
  padj.cutoff   = FDR_cutoff,
  min_count     = min_count,
  min_total     = min_total,
  correction_mode = correction_mode,
  batch_var       = batch_var,
  ruv_k           = ruv_k,
  ruv_controls    = ruv_controls,
  top_n         = top_genes
)

# Save workspace
save(lnccounts, metadata, group, samples, dge_res, file = "DGE_lncRNAs.RData")


r_version    <- paste(R.version\$major, R.version\$minor, sep=".")
edger_version <- as.character(packageVersion("edgeR"))
tidyverse_version <- as.character(packageVersion("tidyverse"))
tibble_version <- as.character(packageVersion("tibble"))
EnhancedVolcano_version <- as.character(packageVersion("EnhancedVolcano"))
RColorBrewer_version <- as.character(packageVersion("RColorBrewer"))
gplots_version <- as.character(packageVersion("gplots"))
openxlsx_version <- as.character(packageVersion("openxlsx"))

versions_text <- paste0(
  "\"", "${task.process}", "\":\n",
  "  r-base: \"", r_version, "\"\n",
  "  edger: \"", edger_version, "\"\n",
  "  tidyverse: \"", tidyverse_version, "\"\n",
  "  tibble: \"", tibble_version, "\"\n",
  "  EnhancedVolcano: \"", EnhancedVolcano_version, "\"\n",
  "  RColorBrewer: \"", RColorBrewer_version, "\"\n",
  "  gplots: \"",gplots_version, "\"\n",
  "  openxlsx: \"", openxlsx_version, "\"\n"

)
writeLines(versions_text, "versions.yml")
