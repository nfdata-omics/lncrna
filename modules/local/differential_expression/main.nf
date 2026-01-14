process DIFFERENTIAL_EXPRESSION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    tuple val(meta), path(count_matrix)
    path gene_info
    path design_file
    val de_method

    output:
    tuple val(meta), path("*.de_results.tsv")       , emit: results
    tuple val(meta), path("*.normalized_counts.tsv"), emit: normalized
    tuple val(meta), path("*.dds.rds")              , emit: dds, optional: true
    path "*.ma_plot.png"                            , emit: ma_plot
    path "*.volcano_plot.png"                       , emit: volcano
    path "*.pca_plot.png"                           , emit: pca
    path "*.heatmap.png"                            , emit: heatmap
    path "*.dispersion_plot.png"                    , emit: dispersion, optional: true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def method = de_method ?: 'DESeq2'
    def padj_cutoff = task.ext.padj_cutoff ?: '0.05'
    def lfc_cutoff = task.ext.lfc_cutoff ?: '1'
    """
    differential_expression.R \\
        --counts $count_matrix \\
        --gene_info $gene_info \\
        --design $design_file \\
        --method $method \\
        --output $prefix \\
        --cores $task.cpus \\
        --padj_cutoff $padj_cutoff \\
        --lfc_cutoff $lfc_cutoff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        deseq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))" 2>/dev/null || echo "not installed")
        edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))" 2>/dev/null || echo "not installed")
        limma: \$(Rscript -e "cat(as.character(packageVersion('limma')))" 2>/dev/null || echo "not installed")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.de_results.tsv
    touch ${prefix}.normalized_counts.tsv
    touch ${prefix}.dds.rds
    touch ${prefix}.ma_plot.png
    touch ${prefix}.volcano_plot.png
    touch ${prefix}.pca_plot.png
    touch ${prefix}.heatmap.png
    touch ${prefix}.dispersion_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: "4.3.0"
        deseq2: "1.40.0"
    END_VERSIONS
    """
}
