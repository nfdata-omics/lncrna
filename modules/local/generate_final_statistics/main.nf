process GENERATE_FINAL_STATISTICS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' :
        'quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882e231a7b5488b5485:47ef981087c59f79fdbcab4d9d7316e9ac2e688d-0' }"

    input:
    tuple val(meta), path(final_gtf)
    tuple val(meta2), path(lncrna_gtf)
    path cpat_lncrna_results
    path cpat_coding_results
    path classification

    output:
    path "*.final_report.html"          , emit: final_stats_report
    path "*.statistics_summary.tsv"     , emit: final_stats_summary
    path "*.cpat_comparison.png"        , emit: cpat_plot
    path "*.classification_barplot.png" , emit: cpat_classification_plot
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_final_statistics.py \\
        --final_gtf $final_gtf \\
        --lncrna_gtf $lncrna_gtf \\
        --cpat_lncrna $cpat_lncrna_results \\
        --cpat_coding $cpat_coding_results \\
        --classification $classification \\
        --prefix $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """
}
