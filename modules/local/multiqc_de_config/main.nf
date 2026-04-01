/*
 * MULTIQC_DE_CONFIG
 * -----------------
 * Generates a multiqc_de.yml dynamically from the contrast string,
 * so each contrast gets its own section in the MultiQC report.
 *
 * Search patterns use filename-based matching (fn / fn_re) instead of
 * full paths, since Nextflow flattens directory structure when staging
 * files into the MULTIQC workdir.
 *
 * Input : contrast string  (e.g. "CaseA/CaseB:control/control")
 * Output: multiqc_de.yml   (passed as custom config to MULTIQC)
 */

process MULTIQC_DE_CONFIG {
    tag "de_config"
    label 'process_single'

    executor 'local'

    input:
    val contrast

    output:
    path "multiqc_de.yml", emit: config

    exec:
    // ----------------------------------------------------------
    // 1. Parse contrast string → list of "case_vs_control" names
    //    Mirrors the R script logic exactly:
    //      - tolower()
    //      - spaces → underscore
    //      - paste0(test, "-", ctrl) → gsub("-","_vs_", ...)
    // ----------------------------------------------------------
    def parts    = contrast.toLowerCase().split(":")
    def tests    = parts[0].split("/").collect { it.trim().replaceAll(" ", "_") }
    def controls = parts[1].split("/").collect { it.trim().replaceAll(" ", "_") }

    def contrastNames = tests.withIndex().collect { test, i ->
        def ctrl = (i < controls.size()) ? controls[i] : controls[-1]
        "${test}_vs_${ctrl}"
    }

    // ----------------------------------------------------------
    // 2. Build YAML
    // ----------------------------------------------------------
    def yaml = new StringBuilder()

    // ── Section order ───────────────────────────────────────────
    yaml << """\
report_section_order:
  "de-global-plots":
    order: -900
  "de-summary-table":
    order: -901
  "de-normalized-cpm":
    order: -902
"""
    contrastNames.eachWithIndex { name, idx ->
        yaml << "  \"de-contrast-${name}\":\n"
        yaml << "    order: ${-910 - idx}\n"
    }

    yaml << "\nexport_plots: true\ndisable_version_detection: true\n"

    // ── custom_data ─────────────────────────────────────────────
    yaml << "\ncustom_data:\n"

    // 2a. Global QC plots — matched by exact filename
    yaml << """\

  de_global_plots:
    id: "de-global-plots"
    section_name: "DE: Global QC Plots"
    description: >
      Global quality-control plots generated before per-contrast testing.
      MDS and PCA show sample-level clustering; the heatmap shows the
      top 50 most variable lncRNAs across all groups.
    plot_type: "image"
    file_format: "png"

"""

    // 2b. DGE summary table
    yaml << """\
  de_summary_table:
    id: "de-summary-table"
    section_name: "DE: DGE Summary"
    description: >
      Number of differentially expressed lncRNAs per contrast (QLF method).
      Values are reported at the LFC and FDR thresholds used in the analysis.
    plot_type: "table"
    file_format: "tsv"
    headers:
      Contrast:
        title: "Contrast"
        description: "Comparison (case vs control)"
      Down:
        title: "Down-regulated"
        description: "Genes with logFC <= -LFC and FDR <= threshold"
        format: "{:,.0f}"
        scale: "Blues"
      NotSig:
        title: "Not significant"
        description: "Genes not passing thresholds"
        format: "{:,.0f}"
      Up:
        title: "Up-regulated"
        description: "Genes with logFC >= LFC and FDR <= threshold"
        format: "{:,.0f}"
        scale: "Reds"
      LFC:
        title: "LFC threshold"
        format: "{:.2f}"
      FDR:
        title: "FDR threshold"
        format: "{:.3f}"

"""

    // 2c. Normalized CPM table
    yaml << """\
  de_normalized_cpm:
    id: "de-normalized-cpm"
    section_name: "DE: Normalized CPM"
    description: >
      TMM-normalized counts per million (CPM) across all retained samples.
    plot_type: "table"
    file_format: "csv"

"""

    // 2d. One image section per contrast
    contrastNames.each { name ->
        def title = name.replaceAll("_", " ").split(" ").collect { it.capitalize() }.join(" ")
        yaml << """\
  de_contrast_${name}:
    id: "de-contrast-${name}"
    section_name: "DE: ${title}"
    description: >
      Volcano plot and heatmap for contrast <strong>${name}</strong>.
    plot_type: "image"
    file_format: "png"

"""
    }

    // ── sp (search patterns) ────────────────────────────────────
    // NOTE: fn / fn_re match against filename only (not full path)
    // because Nextflow flattens directory structure when staging
    // files into MULTIQC's workdir.
    yaml << "sp:\n"

    // Global plots — fn_re grouping all three (fn only accepts a single string in MultiQC v1.33+)
    yaml << """\

  de_global_plots:
    fn_re: ^(MDS|PCA|heatmap_global)\\.png\$

"""

    // Summary TSV — filename contains contrast name
    yaml << """\
  de_summary_table:
    fn_re: ^DGE_summary_.+\\.tsv\$

"""

    // Normalized CPM — filename contains contrast name
    yaml << """\
  de_normalized_cpm:
    fn_re: ^normalized_CPM_.+\\.csv\$

"""

    // Per-contrast plots — volcano and heatmap filenames contain contrast name
    // volcanoplot_<contrast>.png  and  heatmap_<contrast>.png
    contrastNames.each { name ->
        yaml << """\
  de_contrast_${name}:
    fn_re: ^(volcanoplot|heatmap)_${name}\\.png\$

"""
    }

    // ----------------------------------------------------------
    // 3. Write file
    // ----------------------------------------------------------
    def outFile = task.workDir.resolve("multiqc_de.yml")
    outFile.text = yaml.toString()
}
