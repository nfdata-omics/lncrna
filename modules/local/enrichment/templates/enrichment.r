#!/usr/bin/env Rscript

# -------------------------------
# Guilt-by-association enrichment analysis for lncRNAs
# Based on: miRNA targets enrichment script (Karla Ruiz - IU2 OMICS HT)
#
# Strategy:
#   1. Take DE lncRNAs (UP and DOWN separately) from DGE results
#   2. Compute Pearson correlation between each DE lncRNA and all protein-coding
#      mRNAs using log-CPM normalized expression
#   3. Select highly correlated mRNAs (|r| > threshold, p < 0.05)
#   4. Run GO (BP, CC, MF) + KEGG enrichment on those co-expressed mRNAs
#
# This approach is preferred over direct lncRNA enrichment because GO/KEGG
# annotations are built on protein-coding genes. lncRNAs are functionally
# inferred through their co-expressed coding neighbors.
# -------------------------------

suppressMessages(library(ggplot2))
suppressMessages(library(clusterProfiler))
suppressMessages(library(AnnotationHub))
suppressMessages(library(enrichplot))
suppressMessages(library(ggarchery))
suppressMessages(library(openxlsx))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(edgeR))


#----- Input arguments -----
option_list <- list(
  make_option(c("--cor_method"),    action = "store", type = "character", default = "pearson",
              help = "Correlation method: pearson or spearman (default: pearson)"),
  make_option(c("--cor_cutoff"),    action = "store", type = "double",    default = 0.7,
              help = "Absolute correlation threshold (default: 0.7)"),
  make_option(c("--pval_cutoff"),   action = "store", type = "double",    default = 0.05,
              help = "Correlation p-value cutoff (default: 0.05)"),
  make_option(c("--min_genes"),     action = "store", type = "integer",   default = 10,
              help = "Minimum number of co-expressed mRNAs required to run enrichment (default: 10)")
)

parser <- OptionParser(
  usage = "%prog [options] rdata mrna_counts metadata organism",
  option_list = option_list,
  prog = "gba_enrichment_lncrna",
  description = "
-----------------------------------------------------------
  Description: Guilt-by-association enrichment for lncRNAs.
               Finds protein-coding mRNAs co-expressed with
               DE lncRNAs and performs GO + KEGG enrichment.
  Required positional arguments:
     1. DGE_lncRNAs.RData : RData from dge_edger_lncrna.R
     2. mrna_counts.tsv   : Raw counts of protein-coding genes
                            (same format and samples as lncRNA matrix)
     3. metadata.csv      : Same metadata used in DGE analysis
     4. organism          : Organism code (e.g. cfa, hsa, mmu, rno)
-----------------------------------------------------------"
)

arguments <- parse_args(parser, args = commandArgs(trailingOnly = TRUE), positional_arguments = TRUE)
opt <- arguments$options

if (length(arguments$args) < 4) {
  stop("Error: Four positional arguments required: rdata, mrna_counts, metadata, organism.")
}

rdata_file     <- arguments$args[1]
mrna_file      <- arguments$args[2]
metadata_file  <- arguments$args[3]
organism       <- arguments$args[4]

print(paste("RData file      :", rdata_file))
print(paste("mRNA counts     :", mrna_file))
print(paste("Metadata        :", metadata_file))
print(paste("Organism        :", organism))
print(paste("Correlation     :", opt$cor_method, "| cutoff:", opt$cor_cutoff))
print(paste("P-value cutoff  :", opt$pval_cutoff))

load(rdata_file)
# Expects: dge_res (list of contrasts from perform_dge_lncrna)
# Each contrast must have: dge_res[[contrast]]$res_DGE_lab
#                          dge_res[[contrast]]$logcounts  (lncRNA log-CPM)


# ============================================================
#  FUNCTIONS
# ============================================================

# ---- 1. Load organism DB from AnnotationHub ----
species_codes <- c("mmu", "hsa", "rno", "mml", "dre", "dme",
                   "ath", "cfa", "bta", "cel", "gga", "ptr", "xtr")
orgdb_names   <- c(
  "org.Mmu.eg.db.sqlite",
  "org.Hs.eg.db.sqlite",
  "org.Rn.eg.db.sqlite",
  "org.Mm.eg.db.sqlite",
  "org.Dr.eg.db.sqlite",
  "org.Dm.eg.db.sqlite",
  "org.At.tair.db.sqlite",
  "org.Cf.eg.db.sqlite",
  "org.Bt.eg.db.sqlite",
  "org.Ce.eg.db.sqlite",
  "org.Gg.eg.db.sqlite",
  "org.Pt.eg.db.sqlite",
  "org.Xl.eg.db.sqlite"
)

orgdb_ann <- function(organism) {
  hub   <- AnnotationHub()
  index <- match(organism, species_codes)
  if (is.na(index)) stop(paste("Organism not found:", organism))
  db_name      <- orgdb_names[index]
  query_result <- query(hub, db_name)
  if (length(query_result) == 0) stop("Database not found in AnnotationHub.")
  return(query_result[[length(query_result)]])
}

orgdb <- orgdb_ann(organism)
dir.create("./enrichment", showWarnings = FALSE, recursive = TRUE)
write.csv(metadata(orgdb), "./enrichment/orgdb_metadata.csv", row.names = FALSE)


# ---- 2. Read and normalize mRNA count matrix ----
read_mrna_counts <- function(mrna_file, metadata_file) {

  raw <- read.table(mrna_file,
                    header       = TRUE,
                    sep          = "\t",
                    row.names    = 1,
                    check.names  = FALSE,
                    quote        = "",
                    comment.char = "")

  raw <- round(as.matrix(raw))
  storage.mode(raw) <- "integer"

  metadata <- read.delim(metadata_file, sep = ",", row.names = 1, check.names = FALSE)

  # Keep only samples present in both
  common_samples <- intersect(rownames(metadata), colnames(raw))
  if (length(common_samples) == 0) {
    stop("No common samples between mRNA count matrix and metadata.")
  }
  if (length(common_samples) < ncol(raw)) {
    warning(paste("Keeping only", length(common_samples), "common samples between mRNA matrix and metadata."))
  }

  raw      <- raw[, common_samples]
  metadata <- metadata[common_samples, , drop = FALSE]

  cat("mRNA matrix loaded:", nrow(raw), "genes x", ncol(raw), "samples\n")
  return(list(counts = raw, metadata = metadata))
}


# ---- 3. Compute log-CPM for mRNA matrix ----
normalize_mrna <- function(mrna_counts) {
  # Basic TMM normalization + log-CPM, no filtering
  # (all protein-coding genes are kept as correlation universe)
  dge_mrna <- DGEList(counts = mrna_counts)
  keep      <- filterByExpr(dge_mrna)
  dge_mrna  <- dge_mrna[keep, , keep.lib.sizes = FALSE]
  dge_mrna  <- calcNormFactors(dge_mrna, method = "TMM")
  logcpm    <- cpm(dge_mrna, log = TRUE)
  cat("mRNA genes after filtering:", nrow(logcpm), "\n")
  return(logcpm)
}


# ---- 4. Guilt-by-association: correlate lncRNAs vs mRNAs ----
# Returns a data frame of mRNA gene IDs correlated with the input lncRNA set
find_coexpressed_mrnas <- function(lnc_ids, lnc_logcpm, mrna_logcpm,
                                   cor_method, cor_cutoff, pval_cutoff) {

  # Keep only lncRNAs present in logcpm matrix
  lnc_ids <- lnc_ids[lnc_ids %in% rownames(lnc_logcpm)]
  if (length(lnc_ids) == 0) {
    warning("None of the DE lncRNAs found in logcpm matrix.")
    return(character(0))
  }

  # Align samples between lncRNA and mRNA matrices
  common_samples <- intersect(colnames(lnc_logcpm), colnames(mrna_logcpm))
  if (length(common_samples) < 3) {
    warning("Less than 3 common samples between lncRNA and mRNA matrices.")
    return(character(0))
  }

  lnc_mat  <- lnc_logcpm[lnc_ids,  common_samples, drop = FALSE]
  mrna_mat <- mrna_logcpm[,         common_samples, drop = FALSE]

  cat("  Computing correlations:", nrow(lnc_mat), "lncRNAs x",
      nrow(mrna_mat), "mRNAs across", length(common_samples), "samples\n")

  # For each lncRNA, correlate against all mRNAs
  # Collect mRNA IDs that pass thresholds for at least one lncRNA
  coexp_mrnas <- c()

  for (lnc in rownames(lnc_mat)) {
    lnc_vec <- as.numeric(lnc_mat[lnc, ])

    # Compute correlation + p-value for all mRNAs at once
    cor_results <- apply(mrna_mat, 1, function(mrna_vec) {
      ct <- tryCatch(
        cor.test(lnc_vec, mrna_vec, method = cor_method),
        error = function(e) NULL
      )
      if (is.null(ct)) return(c(r = NA, p = NA))
      return(c(r = ct$estimate, p = ct$p.value))
    })

    cor_df <- as.data.frame(t(cor_results))
    colnames(cor_df) <- c("r", "p")
    cor_df$gene <- rownames(cor_df)

    # Filter by thresholds
    sig <- cor_df[!is.na(cor_df$r) &
                    abs(cor_df$r) >= cor_cutoff &
                    cor_df$p     <  pval_cutoff, ]

    coexp_mrnas <- union(coexp_mrnas, sig$gene)
  }

  cat("  Co-expressed mRNAs found:", length(coexp_mrnas), "\n")
  return(coexp_mrnas)
}


# ---- 5. Convert Ensembl IDs to SYMBOL + ENTREZ ----
convert_ids <- function(ensembl_ids, orgdb) {
  tryCatch({
    mapped <- AnnotationDbi::select(orgdb,
                                    keys    = ensembl_ids,
                                    columns = c("SYMBOL", "ENTREZID"),
                                    keytype = "ENSEMBL")
    mapped <- mapped[!is.na(mapped$SYMBOL) & !is.na(mapped$ENTREZID), ]
    mapped <- mapped[!duplicated(mapped$ENSEMBL), ]
    return(mapped)
  }, error = function(e) {
    warning("ID conversion failed: ", e$message)
    return(data.frame(ENSEMBL  = character(),
                      SYMBOL   = character(),
                      ENTREZID = character()))
  })
}


# ---- 6. GO enrichment ----
perform_go_enrichment <- function(gene_symbols, orgdb, ontology,
                                  label, enrichment_dir) {
  if (length(gene_symbols) == 0) {
    message(paste("No genes to test for GO", ontology, label)); return(NULL)
  }

  ego <- tryCatch(
    enrichGO(gene          = gene_symbols,
             keyType       = "SYMBOL",
             OrgDb         = orgdb,
             ont           = ontology,
             pAdjustMethod = "bonferroni",
             pvalueCutoff  = 0.1,
             readable      = TRUE),
    error = function(e) { warning("GO failed: ", e$message); return(NULL) }
  )

  if (is.null(ego) || nrow(ego@result) == 0) {
    message(paste("No significant GO", ontology, "results for", label))
    return(NULL)
  }

  ego2      <- pairwise_termsim(ego)
  plots_dir <- file.path(enrichment_dir, "plots")
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  write.xlsx(ego@result,
             file = file.path(enrichment_dir, paste0("GO_", ontology, "_", label, ".xlsx")))

  plots <- list(
    barplot  = ego %>% mutate(qscore = -log10(p.adjust)) %>%
                 barplot(x = "qscore") + ggtitle(paste("GO", ontology, label)),
    dotplot  = dotplot(ego, showCategory = 20) + ggtitle(paste("GO", ontology, label)),
    treeplot = treeplot(ego2, hclust_method = "average") + ggtitle(paste("GO", ontology, label)),
    emapplot = emapplot(ego2) + ggtitle(paste("GO", ontology, label))
  )

  for (plot_name in names(plots)) {
    ggsave(file.path(plots_dir, paste0(ontology, "_", plot_name, "_", label, ".png")),
           plots[[plot_name]], width = 8, height = 6, dpi = 300)
  }

  return(list(ego = ego, ego2 = ego2, plots = plots))
}


# ---- 7. KEGG enrichment ----
perform_kegg_enrichment <- function(entrez_ids, organism, label, enrichment_dir) {
  if (length(entrez_ids) == 0) {
    message(paste("No genes to test for KEGG", label)); return(NULL)
  }

  kegg <- tryCatch(
    enrichKEGG(gene         = entrez_ids,
               organism     = organism,
               pvalueCutoff = 0.1),
    error = function(e) { warning("KEGG failed: ", e$message); return(NULL) }
  )

  if (is.null(kegg) || nrow(kegg@result) == 0) {
    message(paste("No significant KEGG results for", label)); return(NULL)
  }

  write.xlsx(kegg@result,
             file = file.path(enrichment_dir, paste0("KEGG_", label, ".xlsx")))

  ggsave(file.path(enrichment_dir, paste0("KEGG_dotplot_", label, ".png")),
         dotplot(kegg) + ggtitle(paste("KEGG", label)), height = 7, width = 7, dpi = 600)
  ggsave(file.path(enrichment_dir, paste0("KEGG_barplot_", label, ".png")),
         barplot(kegg) + ggtitle(paste("KEGG", label)), height = 7, width = 7, dpi = 600)

  return(kegg)
}


# ============================================================
#  MAIN EXECUTION
# ============================================================

# Load and normalize mRNA counts
mrna_data   <- read_mrna_counts(mrna_file, metadata_file)
mrna_logcpm <- normalize_mrna(mrna_data$counts)

ontologies         <- c("BP", "CC", "MF")
results_enrichment <- list()

for (contrast_name in names(dge_res)) {
  print(paste("=== Processing contrast:", contrast_name, "==="))

  res_df     <- dge_res[[contrast_name]]$res_DGE_lab
  lnc_logcpm <- dge_res[[contrast_name]]$logcounts

  if (is.null(lnc_logcpm)) {
    warning(paste("logcounts not found for contrast:", contrast_name, "- skipping."))
    next
  }

  # DE lncRNAs split by direction
  lnc_up   <- res_df$lncRNA[res_df$diffexpressed == "UP"]
  lnc_down <- res_df$lncRNA[res_df$diffexpressed == "DOWN"]
  print(paste("  DE lncRNAs - UP:", length(lnc_up), "| DOWN:", length(lnc_down)))

  directions <- list(UP = lnc_up, DOWN = lnc_down)

  ctr_dir <- file.path("./enrichment", contrast_name)
  dir.create(ctr_dir, showWarnings = FALSE, recursive = TRUE)

  results_enrichment[[contrast_name]] <- list()

  for (direction in names(directions)) {
    lnc_ids <- directions[[direction]]

    if (length(lnc_ids) == 0) {
      message(paste("  No", direction, "lncRNAs — skipping.")); next
    }

    print(paste(" Processing direction:", direction))

    # Step 1: find co-expressed mRNAs
    coexp_mrnas <- find_coexpressed_mrnas(
      lnc_ids     = lnc_ids,
      lnc_logcpm  = lnc_logcpm,
      mrna_logcpm = mrna_logcpm,
      cor_method  = opt$cor_method,
      cor_cutoff  = opt$cor_cutoff,
      pval_cutoff = opt$pval_cutoff
    )

    # Export co-expression gene list
    enrich_dir <- file.path(ctr_dir, direction)
    dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

    write.csv(data.frame(ensembl_id = coexp_mrnas),
              file      = file.path(enrich_dir, paste0("coexpressed_mRNAs_", direction, ".csv")),
              row.names = FALSE)

    if (length(coexp_mrnas) < opt$min_genes) {
      warning(paste("Only", length(coexp_mrnas), "co-expressed mRNAs found for",
                    contrast_name, direction, "— below min_genes threshold, skipping enrichment."))
      next
    }

    # Step 2: convert Ensembl → SYMBOL + ENTREZ
    mapped <- convert_ids(coexp_mrnas, orgdb)
    print(paste("  Mapped mRNAs:", nrow(mapped), "/", length(coexp_mrnas)))

    if (nrow(mapped) < opt$min_genes) {
      warning(paste("Too few mRNAs mapped to SYMBOL/ENTREZ for", direction, "— skipping enrichment."))
      next
    }

    # Export mapped gene list
    write.csv(mapped,
              file      = file.path(enrich_dir, paste0("coexpressed_mRNAs_mapped_", direction, ".csv")),
              row.names = FALSE)

    # Step 3: GO enrichment
    go_results <- list()
    for (ont in ontologies) {
      print(paste("   GO", ont, "-", direction))
      go_results[[ont]] <- perform_go_enrichment(
        gene_symbols   = mapped$SYMBOL,
        orgdb          = orgdb,
        ontology       = ont,
        label          = paste0(contrast_name, "_", direction),
        enrichment_dir = enrich_dir
      )
    }

    # Step 4: KEGG enrichment
    print(paste("   KEGG -", direction))
    kegg_result <- perform_kegg_enrichment(
      entrez_ids     = mapped$ENTREZID,
      organism       = organism,
      label          = paste0(contrast_name, "_", direction),
      enrichment_dir = enrich_dir
    )

    results_enrichment[[contrast_name]][[direction]] <- list(
      coexpressed_mrnas = coexp_mrnas,
      mapped            = mapped,
      GO                = go_results,
      KEGG              = kegg_result
    )
  }
}

save(results_enrichment, file = "enrichment_lncRNAs_GBA.RData")
print("Guilt-by-association enrichment analysis complete.")
